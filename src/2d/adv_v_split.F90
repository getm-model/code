#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_v_split - 1D y-advection \label{sec-v-split-adv}
!
! !INTERFACE:
   subroutine adv_v_split(dt,f,Di,adv,V,Do,DV,                          &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxv,dyv,arcd1,                                &
#endif
                          az,au,av,splitfac,scheme,AH,onestep_finalise)
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
   use advection, only: flux
   use advection, only: UPSTREAM,P2,SUPERBEE,MUSCL,P2_PDM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: V,Do,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyv,arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: az,au,av
   integer,intent(in)                         :: scheme
   logical,intent(in),optional                :: onestep_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical            :: use_limiter,use_AH
   integer            :: i,j
   REALTYPE           :: Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_v_split() # ',Ncall
#endif

   use_limiter = (scheme .ne. UPSTREAM)
   use_AH = (AH .gt. _ZERO_)

!$OMP PARALLEL DEFAULT(SHARED)                                  &
!$OMP PARALLEL FIRSTPRIVATE(use_limiter)                        &
!$OMP PARALLEL PRIVATE(i,j,Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd)

! Calculating v-interface fluxes !

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin-1,imax+1
         if (av(i,j).eq.1 .or. (av(i,j).eq.2 .and. (az(i,j).eq.1 .or. az(i,j+1).eq.1))) then
!           Note (KK): exclude advection/diffusion of normal velocity at open bdys
            if (V(i,j) .gt. _ZERO_) then
               fc = f(i,j  )               ! central
               if (scheme .ne. UPSTREAM) then
!                 Note (KK): also fall back to upstream near boundaries
                  use_limiter = (av(i,j-1).eq.1 .or. (av(i,j-1).eq.2 .and. az(i,j).eq.1))
               end if
               if (use_limiter) then
                  cfl = splitfac*V(i,j)/DV(i,j)*dt/DYV
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
                  use_limiter = (av(i,j+1).eq.1 .or. (av(i,j+1).eq.2 .and. az(i,j+1).eq.1))
               end if
               if (use_limiter) then
                  cfl = -splitfac*V(i,j)/DV(i,j)*dt/DYV
                  fu = f(i,j+2)            ! upstream
                  fd = f(i,j  )            ! downstream
                  if (abs(fc-fd) .gt. 1.d-10) then
                     r = (fu-fc)/(fc-fd)   ! slope ratio
                  else
                     r = (fu-fc)*1.d10
                  end if
               end if
            end if
            flux(i,j) = V(i,j)*fc
            if (use_limiter) then
               select case (scheme)
                  case ((P2),(P2_PDM))
                     x = one6th*(_ONE_-_TWO_*cfl)
                     Phi = (_HALF_+x) + (_HALF_-x)*r
                     if (scheme.eq.P2) then
!                       KK-TODO: for r=0 limit might be non-zero?!
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
               flux(i,j) = flux(i,j) + V(i,j)*_HALF_*limit*(_ONE_-cfl)*(fd-fc)
            end if
            if (use_AH) then
!              Horizontal diffusion
               flux(i,j) = flux(i,j) - AH*DV(i,j)*(f(i,j+1)-f(i,j  ))/DYV
            end if
         else if (use_AH .and. av(i,j).eq.2 .and. (az(i,j).eq.2 .or. az(i,j+1).eq.2)) then
!           Note (KK): special handling for advection/diffusion of normal velocity at open bdys
!                      (advection/diffusion of tracers near open bdys already included in former case)
!                      outflow condition implies no advection across open bdy
            if (az(i,j) .eq. 2) then ! eastern open bdy (az(i,j).ne.1.and.az(i,j+1).ne.1.and.av(i,j-1).eq.1)
               flux(i,j) = AH*DV(i,j-1)*(f(i,j  )-f(i,j-1))/DYVJM1 ! not used for j=jmin-1
            else ! western open bdy (az(i,j).ne.1.and.az(i,j+1).ne.1.and.av(i,j+1).eq.1)
               flux(i,j) = AH*DV(i,j+1)*(f(i,j+2)-f(i,j+1))/DYVJP1 ! not used for j=jmax
            end if
         else
            flux(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!  Doing the v-advection step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax+1
         if (az(i,j).eq.1 .or. (use_AH .and. az(i,j).eq.2 .and. (av(i,j-1).eq.1 .or. av(i,j).eq.1))) then
!           Note (KK): exclude advection/diffusion of tracers at open bdy cells
!                      special handling for advection/diffusion of normal velocity at open bdys
!                      vanishing diffusive flux across exterior interface must be explicitely prescribed
!                      diffusive flux across interior interface is isolated at exterior flux
            Dio = Di(i,j)
            if (az(i,j) .eq. 1) then
               Di(i,j) =  Dio - splitfac*dt*( V(i,j  )*DXV           &
                                             -V(i,j-1)*DXVJM1)*ARCD1
               advn = splitfac*( flux(i,j  )*DXV           &
                                -flux(i,j-1)*DXVJM1)*ARCD1
            else if (av(i,j-1) .eq. 1) then ! eastern open bdy
!              Note (KK): interior advection/diffusion already included in first case
               advn = -splitfac*flux(i,j  )*DXVJM1*ARCD1
            else ! western open bdy (av(i,j) .eq. 1)
!              Note (KK): interior advection/diffusion already included in first case
               advn =  splitfac*flux(i,j-1)*DXV*ARCD1
            end if
            adv(i,j) = adv(i,j) + advn
            if (.not. present(onestep_finalise)) then
!              Note (KK): do the splitting step
               f(i,j) = ( Dio*f(i,j) - dt*advn ) / Di(i,j)
            end if
         end if
         if (present(onestep_finalise)) then
!           Note (KK): do not update f in case of onestep_finalise=.false. !!!
            if (onestep_finalise) then
               if (az(i,j).eq.1 .or. (az(i,j).eq.2 .and. (    au(i-1,j  ).eq.1 &
                                                          .or.au(i  ,j  ).eq.1 &
                                                          .or.av(i  ,j-1).eq.1 &
                                                          .or.av(i  ,j  ).eq.1 ))) then
!                 Note (KK): exclude tracer open bdy cells but
!                            include all velocity open bdys
                  f(i,j) = ( Do(i,j)*f(i,j) - dt*adv(i,j) ) / Di(i,j)
               end if
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
