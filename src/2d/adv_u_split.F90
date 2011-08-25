#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_u_split - 1D x-advection \label{sec-u-split-adv}
!
! !INTERFACE:
   subroutine adv_u_split(dt,f,U,DU, &
                          delxu,delyu,area_inv,au,splitfac,method,az,AH)
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
! The selector for the schemes is {\tt method}:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt method = UPSTREAM\_SPLIT}: & first-order upstream (monotone) \\
! {\tt method = P2}: & third-order polynomial (non-monotone) \\
! {\tt method = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! {\tt method = MUSCL}: & second-order TVD (monotone) \\
! {\tt method = Superbee}: & second-order TVD (monotone) \\
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
   use advection, only: Di,Dio,cu
   use advection, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection, only: one6th
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,DU
   REALTYPE,dimension(E2DFIELD),intent(in)    :: delxu,delyu,area_inv
   integer,dimension(E2DFIELD),intent(in)     :: au,az
   integer,intent(in)                         :: method
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,j
   REALTYPE        :: c,x,r,Phi,limit,fu,fc,fd
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_u_split() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c,x,r,Phi,limit,fu,fc,fd)

! Calculating u-interface fluxes !
   select case (method)
      case (UPSTREAM_SPLIT)
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin-1,imax
               if (au(i,j) .gt. 0) then
                 if (U(i,j) .gt. _ZERO_) then
                     cu(i,j)=U(i,j)*f(i,j)
                  else
                     cu(i,j)=U(i,j)*f(i+1,j)
                  end if
               else
                  cu(i,j) = _ZERO_
               end if
            end do
         end do
!$OMP END DO NOWAIT
!$OMP BARRIER
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin-1,imax
               if (au(i,j) .gt. 0) then
                  if (U(i,j) .gt. _ZERO_) then
                     c=U(i,j)/DU(i,j)*dt/delxu(i,j)
                     if (au(i-1,j) .gt. 0) then
                        fu=f(i-1,j)         ! upstream
                     else
                        fu=f(i  ,j)
                     end if
                     fc=f(i  ,j)            ! central
                     fd=f(i+1,j)            ! downstream
                     if (abs(fd-fc) .gt. 1.d-10) then
                        r=(fc-fu)/(fd-fc)     ! slope ratio
                     else
                        r=(fc-fu)*1.d10
                     end if
                  else
                     c=-U(i,j)/DU(i,j)*dt/delxu(i,j)
                     if (au(i+1,j) .gt. 0) then
                        fu=f(i+2,j)         ! upstream
                     else
                        fu=f(i+1,j)
                     end if
                     fc=f(i+1,j)            ! central
                     fd=f(i  ,j)            ! downstream
                     if (abs(fc-fd) .gt. 1.d-10) then
                        r=(fu-fc)/(fc-fd)     ! slope ratio
                     else
                        r=(fu-fc)*1.d10
                     end if
                  end if
                  select case (method)
                     case ((P2),(P2_PDM))
                        x = one6th*(_ONE_-_TWO_*c)
                        Phi=(_HALF_+x)+(_HALF_-x)*r
                        if (method.eq.P2) then
                           limit=Phi
                        else
                           limit=max(_ZERO_,min(Phi,_TWO_/(_ONE_-c), &
                                     _TWO_*r/(c+1.d-10)))
                        end if
                     case (Superbee)
                        limit=max(_ZERO_,min(_ONE_,_TWO_*r),min(r,_TWO_))
                     case (MUSCL)
                        limit=max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                     case default
                        FATAL 'Not so good - do_advection()'
                        stop 'adv_u_split'
                  end select
                  cu(i,j)=U(i,j)*(fc+_HALF_*limit*(_ONE_-c)*(fd-fc))
!Horizontal diffusion
                  if ( AH.gt._ZERO_ .and. az(i,j).gt.0 .and. az(i+1,j).gt.0 ) then
                     cu(i,j) = cu(i,j)-AH*DU(i,j) &
                                  *(f(i+1,j)-f(i,j))/delxu(i,j)
                  end if
               else
                  cu(i,j) = _ZERO_
               end if
            end do
         end do
!$OMP END DO NOWAIT
!$OMP BARRIER
   end select

!$OMP BARRIER
!  Doing the u-advection step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            Dio(i,j)=Di(i,j)
            Di(i,j)=Dio(i,j)                           &
                     -splitfac*dt*(U(i,j)*delyu(i,j)    &
                                  -U(i-1,j)*delyu(i-1,j))*area_inv(i,j)
            f(i,j)=(f(i,j)*Dio(i,j)           &
                     -splitfac*dt*(cu(i,j)*delyu(i,j)    &
                     -cu(i-1,j)*delyu(i-1,j))*area_inv(i,j) &
                     )/Di(i,j)
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_u_split()'
   write(debug,*)
#endif
   return
   end subroutine adv_u_split

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
