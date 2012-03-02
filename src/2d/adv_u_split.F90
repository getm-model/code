#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_u_split - 1D x-advection \label{sec-u-split-adv}
!
! !INTERFACE:
   subroutine adv_u_split(dt,f,Di,adv,U,Do,DU,            &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxu,dyu,arcd1,                  &
#endif
                          splitfac,scheme,AH,             &
                          mask_flux,mask_update,          &
                          nosplit_finalise,mask_finalise)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! Here, the $x$-directional split 1D advection step is executed
! with a number of options for the numerical scheme.
!
! Here, $n$ and $o$ denote values before and after this operation,
! respectively, $n$ denote intermediate values when other
! 1D advection steps come after this and $o$ denotes intermediate
! values when other 1D advection steps came before this.
! Furthermore, when this $x$-directional split step is repeated
! during the total time step (Strang splitting), the time step $\Delta t$
! denotes a fraction of the full time step.
!
! The interfacial fluxes are calculated by means of
! monotone Total Variation Diminishing (TVD), the first-order monotone
! upstream and the (non-monotone)
! unlimited third-order polynomial scheme.
!
! The selector for the schemes is {\tt scheme}:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt scheme = UPSTREAM}: & first-order upstream (monotone) \\
! {\tt scheme = P2}: & third-order polynomial (non-monotone) \\
! {\tt scheme = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! {\tt scheme = MUSCL}: & second-order TVD (monotone) \\
! {\tt scheme = SUPERBEE}: & second-order TVD (monotone) \\
! \end{tabular}
!
! \vspace{0.5cm}
!
! Furthermore, the horizontal diffusion in $x$-direction
! with the constant diffusion
! coefficient {\tt AH} is carried out here by means of a central difference
! second-order scheme.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use advection, only: uflux
   use advection, only: UPSTREAM,P2,SUPERBEE,MUSCL,P2_PDM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  Note (KK): in general dxu, dyu and mask_flux do only have valid data
!             within (_IRANGE_HALO_-1,_JRANGE_HALO_). In some cases the
!             original field extension may even be _IRANGE_HALO_. Then
!             explicit declared array bounds _IRANGE_HALO_-1 require a
!             provision of the corresponding subarray and will cause
!             copying of the non-contiguously data into a temporarily
!             array. Therefore they are declared as pointers here. This
!             however requires, that the provided pointers already carry
!             the correct bounds.
   REALTYPE,intent(in)                             :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)         :: U,Do,DU
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in)      :: dxu,dyu
   REALTYPE,dimension(E2DFIELD),intent(in)         :: arcd1
#endif
   integer,intent(in)                              :: scheme
   logical,dimension(:,:),pointer,intent(in)       :: mask_flux
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
   integer            :: i,j,iadd
   REALTYPE           :: dti,Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_u_split() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif

   if (scheme .eq. UPSTREAM) then
      iadd = 1
   else
      iadd = 0
   end if

   use_limiter = .false.
   use_AH = (AH .gt. _ZERO_)
   dti = splitfac*dt

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j,use_limiter)                             &
!$OMP          PRIVATE(i,Dio,advn,cfl,x,r,Phi,limit,fu,fc,fd)

! Calculating u-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-1-iadd,imax+iadd
         if (mask_flux(i,j)) then
!           Note (KK): exclude x-advection of u across W/E open bdys
            if (U(i,j) .gt. _ZERO_) then
               fc = f(i  ,j)               ! central
               if (scheme .ne. UPSTREAM) then
!                 Note (KK): also fall back to upstream near boundaries
                  use_limiter = mask_flux(i-1,j)
               end if
               if (use_limiter) then
                  cfl = U(i,j)/DU(i,j)*dti/DXU
                  fu = f(i-1,j)            ! upstream
                  fd = f(i+1,j)            ! downstream
                  if (abs(fd-fc) .gt. 1.d-10) then
                     r = (fc-fu)/(fd-fc)   ! slope ratio
                  else
                     r = (fc-fu)*1.d10
                  end if
               end if
            else
               fc = f(i+1,j)               ! central
               if (scheme .ne. UPSTREAM) then
!                 Note (KK): also fall back to upstream near boundaries
                  use_limiter = mask_flux(i+1,j)
               end if
               if (use_limiter) then
                  cfl = -U(i,j)/DU(i,j)*dti/DXU
                  fu = f(i+2,j)            ! upstream
                  fd = f(i  ,j)            ! downstream
                  if (abs(fc-fd) .gt. 1.d-10) then
                     r = (fu-fc)/(fc-fd)   ! slope ratio
                  else
                     r = (fu-fc)*1.d10
                  end if
               end if
            end if
            uflux(i,j) = U(i,j)*fc
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
                     stop 'adv_u_split: invalid scheme'
               end select
               uflux(i,j) = uflux(i,j) + U(i,j)*_HALF_*limit*(_ONE_-cfl)*(fd-fc)
            end if
            if (use_AH) then
!              Horizontal diffusion
               uflux(i,j) = uflux(i,j) - AH*DU(i,j)*(f(i+1,j)-f(i  ,j))/DXU
            end if
         else
            uflux(i,j) = _ZERO_
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-iadd,imax+iadd
         if (mask_update(i,j)) then
!           Note (KK): exclude x-advection of tracer and u across W/E open bdys
            Dio = Di(i,j)
            Di(i,j) =  Dio - dti*( U(i  ,j)*DYU           &
                                  -U(i-1,j)*DYUIM1)*ARCD1
            advn = splitfac*( uflux(i  ,j)*DYU           &
                             -uflux(i-1,j)*DYUIM1)*ARCD1
            adv(i,j) = adv(i,j) + advn
            if (.not. present(nosplit_finalise)) then
!              do the x-advection splitting step
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
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_u_split()'
   write(debug,*)
#endif
   return
   end subroutine adv_u_split
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
