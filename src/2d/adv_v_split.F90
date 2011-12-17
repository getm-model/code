#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_v_split - 1D y-advection \label{sec-v-split-adv}
!
! !INTERFACE:
   subroutine adv_v_split(dt,f,Di,adv,V,Do,DV,            &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxv,dyv,arcd1,                  &
#endif
                          splitfac,scheme,AH,             &
                          mask_flux,mask_update,          &
                          nosplit_finalise,mask_finalise)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! Here, the $y$-directional split 1D advection step is executed
! with a number of options for the numerical scheme. The basic
! advection equation is accompanied by an fractional step
! for the continuity equation.
!
! When this $y$-directional split step is repeated
! during the total time step (Strang splitting), the time step $\Delta t$
! denotes a fraction of the full time step.
!
! The interfacial fluxes $\tilde c^v_{i,j,k}$ are calculated by means of
! monotone and non-monotone schemes which are described in detail in
! {\tt u\_split\_adv}, see section \ref{sec-u-split-adv} on
! page \pageref{sec-u-split-adv}.
!
! Furthermore, the horizontal diffusion in $y$-direction
! with the constant diffusion
! coefficient {\tt AH} is carried out here by means of a central difference
! second-order scheme.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use advection, only: vflux
   use advection, only: UPSTREAM,P2,SUPERBEE,MUSCL,P2_PDM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                             :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)         :: V,Do,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)         :: arcd1
#endif
   integer,intent(in)                              :: scheme
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: mask_flux
   logical,dimension(E2DFIELD),intent(in)          :: mask_update
   logical,intent(in),optional                     :: nosplit_finalise
   logical,dimension(E2DFIELD),intent(in),optional :: mask_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)      :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical            :: use_limiter,use_AH
   integer            :: i,j,jadd
   REALTYPE           :: dti,Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_v_split() # ',Ncall
#endif

   if (scheme .eq. UPSTREAM) then
      jadd = 1
   else
      jadd = 0
   end if

   use_limiter = .false.
   use_AH = (AH .gt. _ZERO_)
   dti = splitfac*dt

!$OMP PARALLEL DEFAULT(SHARED)                                  &
!$OMP          FIRSTPRIVATE(use_limiter)                        &
!$OMP          PRIVATE(i,j,Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd)

! Calculating v-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1-jadd,jmax+jadd
      do i=imin-HALO,imax+HALO
         if (mask_flux(i,j)) then
!           Note (KK): exclude y-advection of v across N/S open bdys
            if (V(i,j) .gt. _ZERO_) then
               fc = f(i,j  )               ! central
               if (scheme .ne. UPSTREAM) then
!                 Note (KK): also fall back to upstream near boundaries
                  use_limiter = mask_flux(i,j-1)
               end if
               if (use_limiter) then
                  cfl = V(i,j)/DV(i,j)*dti/DYV
                  fu = f(i,j-1)            ! upstream
                  fd = f(i,j+1)            ! downstream
                  if (abs(fd-fc) .gt. 1.d-10) then
                     r = (fc-fu)/(fd-fc)   ! slope ratio
                  else
                     r = (fc-fu)*1.d10
                  end if
               end if
            else
               fc = f(i,j+1)               ! central
               if (scheme .ne. UPSTREAM) then
!                 Note (KK): also fall back to upstream near boundaries
                  use_limiter = mask_flux(i,j+1)
               end if
               if (use_limiter) then
                  cfl = -V(i,j)/DV(i,j)*dti/DYV
                  fu = f(i,j+2)            ! upstream
                  fd = f(i,j  )            ! downstream
                  if (abs(fc-fd) .gt. 1.d-10) then
                     r = (fu-fc)/(fc-fd)   ! slope ratio
                  else
                     r = (fu-fc)*1.d10
                  end if
               end if
            end if
            vflux(i,j) = V(i,j)*fc
            if (use_limiter) then
               select case (scheme)
                  case ((P2),(P2_PDM))
                     x = one6th*(_ONE_-_TWO_*cfl)
                     Phi = (_HALF_+x) + (_HALF_-x)*r
                     if (scheme.eq.P2) then
                        limit = Phi
                     else
                        limit = max(_ZERO_,min(Phi,_TWO_/(_ONE_-cfl),_TWO_*r/(cfl+1.d-10)))
                     end if
                  case (SUPERBEE)
                     limit = max(_ZERO_,min(_ONE_,_TWO_*r),min(r,_TWO_))
                  case (MUSCL)
                     limit = max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                  case default
                     stop 'adv_v_split: invalid scheme'
               end select
               vflux(i,j) = vflux(i,j) + V(i,j)*_HALF_*limit*(_ONE_-cfl)*(fd-fc)
            end if
            if (use_AH) then
!              Horizontal diffusion
               vflux(i,j) = vflux(i,j) - AH*DV(i,j)*(f(i,j+1)-f(i,j  ))/DYV
            end if
         else
            vflux(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-jadd,jmax+jadd
      do i=imin-HALO,imax+HALO
         if (mask_update(i,j)) then
!           Note (KK): exclude y-advection of tracer and v across N/S open bdys
            Dio = Di(i,j)
            Di(i,j) =  Dio - dti*( V(i,j  )*DXV           &
                                  -V(i,j-1)*DXVJM1)*ARCD1
            advn = splitfac*( vflux(i,j  )*DXV           &
                             -vflux(i,j-1)*DXVJM1)*ARCD1
            adv(i,j) = adv(i,j) + advn
            if (.not. present(nosplit_finalise)) then
!              do the y-advection splitting step
               f(i,j) = ( Dio*f(i,j) - dt*advn ) / Di(i,j)
            end if
         end if
         if (present(nosplit_finalise)) then
            if (nosplit_finalise .and. mask_finalise(i,j)) then
!              Note (KK): do not modify tracer inside open bdy cells
               f(i,j) = ( Do(i,j)*f(i,j) - dt*adv(i,j) ) / Di(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_v_split()'
   write(debug,*)
#endif
   return
   end subroutine adv_v_split
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
