#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_v_split - 1D y-advection \label{sec-v-split-adv}
!
! !INTERFACE:
   subroutine adv_v_split(dt,f,V,DV, &
                          delxv,delyv,area_inv,av,splitfac,method,az,AH)
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
   use advection, only: Di,Dio,cu
   use advection, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection, only: one6th
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: V,DV
   REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyv,area_inv
   integer,dimension(E2DFIELD),intent(in)     :: av,az
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
   write(debug,*) 'adv_v_split() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c,x,r,Phi,limit,fu,fc,fd)

! Calculating v-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-1,jmax
            do i=imin,imax
               if(av(i,j) .gt. 0) then
                  if (V(i,j) .gt. _ZERO_) then
                     cu(i,j)=V(i,j)*f(i,j)
                  else
                     cu(i,j)=V(i,j)*f(i,j+1)
                  end if
               else
                  cu(i,j) = _ZERO_
               end if
            end do
         end do
!$OMP END DO NOWAIT
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-1,jmax
            do i=imin,imax
               if (av(i,j) .gt. 0) then
                  if (V(i,j) .gt. _ZERO_) then
                     c=V(i,j)/DV(i,j)*dt/delyv(i,j)
                     if (av(i,j-1) .gt. 0) then
                        fu=f(i,j-1)         ! upstream
                     else
                        fu=f(i,j  )
                     end if
                     fc=f(i,j  )            ! central
                     fd=f(i,j+1)            ! downstream
                     if (abs(fd-fc) .gt. 1.d-10) then
                        r=(fc-fu)/(fd-fc)     ! slope ratio
                     else
                        r=(fc-fu)*1.d10
                     end if
                  else
                     c=-V(i,j)/DV(i,j)*dt/delyv(i,j)
                     if (av(i,j+1) .gt. 0) then
                        fu=f(i,j+2)         ! upstream
                     else
                        fu=f(i,j+1)
                     end if
                     fc=f(i,j+1)            ! central
                     fd=f(i,j  )            ! downstream
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
                           limit=max(_ZERO_,min(Phi,_TWO_/(_ONE_-c),_TWO_*r/(c+1.d-10)))
                        end if
                     case (Superbee)
                        limit=max(_ZERO_, min(_ONE_,_TWO_*r), min(r,_TWO_) )
                     case (MUSCL)
                        limit=max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                     case default
                        FATAL 'This is not so good - do_advection()'
                        stop 'adv_v_split'
                  end select
                  cu(i,j)=V(i,j)*(fc+_HALF_*limit*(_ONE_-c)*(fd-fc))
!Horizontal diffusion
                  if ( AH.gt._ZERO_ .and. az(i,j).gt.0 .and. az(i,j+1).gt.0 ) &
                     cu(i,j)=cu(i,j)-AH*DV(i,j)*(f(i,j+1)-f(i,j))/delyv(i,j)
               else
                  cu(i,j) = _ZERO_
               end if
            end do
         end do
!$OMP END DO NOWAIT
   end select

!$OMP BARRIER
!  Doing the v-advection step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j).eq.1) then
            Dio(i,j)=Di(i,j)
            Di(i,j)=Dio(i,j)              &
                     -splitfac*dt*(V(i,j)*delxv(i,j)    &
                     -V(i,j-1)*delxv(i,j-1))*area_inv(i,j)
            f(i,j)=(f(i,j)*Dio(i,j)        &
                    -splitfac*dt*(cu(i,j)*delxv(i,j)    &
                    -cu(i,j-1)*delxv(i,j-1))*area_inv(i,j))/Di(i,j)
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_v_split()'
   write(debug,*)
#endif
   return
   end subroutine adv_v_split

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
