#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_u_split - 1D x-advection \label{sec-u-split-adv}
!
! !INTERFACE:
   subroutine adv_u_split(dt,f,Di,adv,U,Do,DU,                           &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxu,dyu,arcd1,                                 &
#endif
                          au,splitfac,adv_scheme,az,AH,onestep_finalise)
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
! The selector for the schemes is {\tt adv_scheme}:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt adv_scheme = UPSTREAM}: & first-order upstream (monotone) \\
! {\tt adv_scheme = P2}: & third-order polynomial (non-monotone) \\
! {\tt adv_scheme = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! {\tt adv_scheme = MUSCL}: & second-order TVD (monotone) \\
! {\tt adv_scheme = SUPERBEE}: & second-order TVD (monotone) \\
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
   use advection, only: cu
   use advection, only: UPSTREAM,P2,SUPERBEE,MUSCL,P2_PDM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,Do,DU
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(E2DFIELD),intent(in)    :: dxu,dyu,arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: au,az
   integer,intent(in)                         :: adv_scheme
   logical,intent(in),optional                :: onestep_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer            :: i,j
   REALTYPE           :: Dio,advn,c,x,r,Phi,limit,fu,fc,fd
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_u_split() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,Dio,advn,c,x,r,Phi,limit,fu,fc,fd)

! Calculating u-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (au(i,j).eq.1 .or. (au(i,j).eq.2 .and. (az(i,j).eq.1 .or. az(i+1,j).eq.1))) then
!           Note (KK): exclude advection/diffusion of normal velocity at open bdys
            if (U(i,j) .gt. _ZERO_) then
               fc=f(i  ,j)            ! central
               if (adv_scheme .ne. UPSTREAM) then
                  c=U(i,j)/DU(i,j)*dt/DXU
                  if (au(i-1,j) .gt. 0) then
                     fu=f(i-1,j)         ! upstream
                  else
                     fu=f(i  ,j)
                  end if
                  fd=f(i+1,j)            ! downstream
                  if (abs(fd-fc) .gt. 1.d-10) then
                     r=(fc-fu)/(fd-fc)     ! slope ratio
                  else
                     r=(fc-fu)*1.d10
                  end if
               end if
            else
               fc=f(i+1,j)            ! central
               if (adv_scheme .ne. UPSTREAM) then
                  c=-U(i,j)/DU(i,j)*dt/DXU
                  if (au(i+1,j) .gt. 0) then
                     fu=f(i+2,j)         ! upstream
                  else
                     fu=f(i+1,j)
                  end if
                  fd=f(i  ,j)            ! downstream
                  if (abs(fc-fd) .gt. 1.d-10) then
                     r=(fu-fc)/(fc-fd)     ! slope ratio
                  else
                     r=(fu-fc)*1.d10
                  end if
               end if
            end if
            cu(i,j) = U(i,j)*fc
            if (adv_scheme .ne. UPSTREAM) then
               select case (adv_scheme)
                  case ((P2),(P2_PDM))
                     x = one6th*(_ONE_-_TWO_*c)
                     Phi=(_HALF_+x)+(_HALF_-x)*r
                     if (adv_scheme.eq.P2) then
                        limit=Phi
                     else
                        limit=max(_ZERO_,min(Phi,_TWO_/(_ONE_-c),_TWO_*r/(c+1.d-10)))
                     end if
                  case (SUPERBEE)
                     limit=max(_ZERO_,min(_ONE_,_TWO_*r),min(r,_TWO_))
                  case (MUSCL)
                     limit=max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                  case default
                     FATAL 'Invalid adv_scheme'
                     stop
               end select
               cu(i,j) = cu(i,j) + U(i,j)*_HALF_*limit*(_ONE_-c)*(fd-fc)
            end if
            if (AH .gt. _ZERO_) then
!              Horizontal diffusion
               cu(i,j) = cu(i,j) - AH*DU(i,j)*(f(i+1,j)-f(i  ,j))/DXU
            end if
         else if (AH.gt._ZERO_ .and. au(i,j).eq.2 .and. (az(i,j).eq.2 .or. az(i+1,j).eq.2)) then
!           Note (KK): special handling for advection/diffusion of normal velocity at open bdys
!                      (advection/diffusion of tracers at open bdys already included in former case)
!                      outflow condition implies no advection across open bdy
            if (az(i,j) .eq. 2) then ! eastern open bdy (az(i,j).ne.1.and.az(i+1,j).ne.1.and.au(i-1,j).eq.1)
               cu(i,j) = AH*DU(i-1,j)*(f(i  ,j)-f(i-1,j))/DXUIM1
            else ! western open bdy (az(i,j).ne.1.and.az(i+1,j).ne.1.and.au(i+1,j).eq.1)
               cu(i,j) = AH*DU(i+1,j)*(f(i+2,j)-f(i+1,j))/DXUIP1
            end if
         else
            cu(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!  Doing the u-advection step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j).eq.1 .or. (AH.gt._ZERO_ .and. az(i,j).eq.2 .and. (au(i-1,j).eq.1 .or. au(i,j).eq.1))) then
!           Note (KK): exclude advection/diffusion of tracers at open bdy cells
!                      special handling for advection/diffusion of normal velocity at open bdys
!                      vanishing diffusive flux across exterior interface must be explicitely prescribed
!                      diffusive flux across interior interface is isolated at exterior cu
            Dio = Di(i,j)
            if (az(i,j) .eq. 1) then
               Di(i,j) =  Dio - splitfac*dt*( U(i  ,j)*DYU           &
                                             -U(i-1,j)*DYUIM1)*ARCD1
               advn = splitfac*( cu(i  ,j)*DYU           &
                                -cu(i-1,j)*DYUIM1)*ARCD1
            else if (au(i-1,j) .eq. 1) then ! eastern open bdy
!              Note (KK): interior advection/diffusion already included in first case
               advn = -splitfac*cu(i  ,j)*DYUIM1*ARCD1
            else ! western open bdy (au(i,j) .eq. 1)
!              Note (KK): interior advection/diffusion already included in first case
               advn =  splitfac*cu(i-1,j)*DYU*ARCD1
            end if
            adv(i,j) = adv(i,j) + advn
            if (present(onestep_finalise)) then
!              Note (KK): do not update f in case of onestep_finalise=.false. !!!
               if (onestep_finalise) then
                  f(i,j) = ( Do(i,j)*f(i,j) - dt*adv(i,j) ) / Di(i,j)
               end if
            else
               f(i,j) = ( Dio*f(i,j) - dt*advn ) / Di(i,j)
            end if
         end if
      end do
   end do
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
