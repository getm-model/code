#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  v_split_adv - 1D y-advection \label{sec-v-split-adv}
!
! !INTERFACE:
   subroutine v_split_adv(dt,f,vv,hvn, &
                          delxv,delyv,area_inv,av,splitfac,method,az)
! !DESCRIPTION:
!
! Here, the $y$-directional split 1D advection step is executed
! with a number of options for the numerical scheme. The basic
! advection equation is accompanied by an fractional step
! for the continuity equation and both equations look as follows:
!
! \begin{equation}\label{adv_v_step}
! h^n_{i,j,k} c^n_{i,j,k} =
! h^o_{i,j,k} c^o_{i,j,k}
! \displaystyle
! - \Delta t
! \frac{
! q_{i,j,k}\tilde c^v_{i,j,k}\Delta y^v_{i,j}-
! q_{i,j-1,k}\tilde c^v_{i,j-1,k}\Delta y^v_{i,j-1}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}},
! \end{equation}
!
! with the 1D continuity equation
! \begin{equation}\label{adv_v_step_h}
! h^n_{i,j,k} =
! h^o_{i,j,k}
! \displaystyle
! - \Delta t
! \frac{
! q_{i,j,k}\Delta y^v_{i,j}-
! q_{i,j-1,k}\Delta y^v_{i,j-1}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}.
! \end{equation}
!
! Here, $n$ and $o$ denote values before and after this operation,
! respectively, $n$ denote intermediate values when other
! 1D advection steps come after this and $o$ denotes intermediate
! values when other 1D advection steps came before this.
! Furthermore, when this $y$-directional split step is repeated
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
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection_3d, only: hi,hio,cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: one6th
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: vv(I3DFIELD),hvn(I3DFIELD)
   REALTYPE, intent(in) :: delxv(I2DFIELD),delyv(I2DFIELD)
   REALTYPE, intent(in) :: area_inv(I2DFIELD),dt
   integer, intent(in)  :: av(E2DFIELD),az(E2DFIELD)
   REALTYPE, intent(in) :: splitfac
   integer, intent(in)  :: method
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,j,k
   REALTYPE        :: c,x,r,Phi,limit,fu,fc,fd
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'v_split_adv() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,c,x,r,Phi,limit,fu,fc,fd)

! Calculating v-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-1,jmax
               do i=imin,imax
                  if(av(i,j) .gt. 0) then
                     if (vv(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=vv(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=vv(i,j,k)*f(i,j+1,k)
                     end if
                  else
                     cu(i,j,k) = _ZERO_
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-1,jmax
               do i=imin,imax
                  if (av(i,j) .gt. 0) then
                     if (vv(i,j,k) .gt. _ZERO_) then
                        c=vv(i,j,k)/hvn(i,j,k)*dt/delyv(i,j)
                        if (av(i,j-1) .gt. 0) then
                           fu=f(i,j-1,k)         ! upstream
                        else
                           fu=f(i,j  ,k)
                        end if
                        fc=f(i,j  ,k)            ! central
                        fd=f(i,j+1,k)            ! downstream
                        if (abs(fd-fc) .gt. 1.d-10) then
                           r=(fc-fu)/(fd-fc)     ! slope ratio
                        else
                           r=(fc-fu)*1.d10
                        end if
                     else
                        c=-vv(i,j,k)/hvn(i,j,k)*dt/delyv(i,j)
                        if (av(i,j+1) .gt. 0) then
                           fu=f(i,j+2,k)         ! upstream
                        else
                           fu=f(i,j+1,k)
                        end if
                        fc=f(i,j+1,k)            ! central
                        fd=f(i,j  ,k)            ! downstream
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
                           FATAL 'This is not so good - do_advection_3d()'
                           stop 'v_split_adv'
                     end select
                     cu(i,j,k)=vv(i,j,k)*(fc+_HALF_*limit*(_ONE_-c)*(fd-fc))
                  else
                     cu(i,j,k) = _ZERO_
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do
   end select

!$OMP BARRIER
   do k=1,kmax   ! Doing the v-advection step
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j).eq.1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)              &
                        -splitfac*dt*(vv(i,j,k)*delxv(i,j)    &
                        -vv(i,j-1,k)*delxv(i,j-1))*area_inv(i,j)
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)        &
                       -splitfac*dt*(cu(i,j,k)*delxv(i,j)    &
                       -cu(i,j-1,k)*delxv(i,j-1))*area_inv(i,j))/hi(i,j,k)
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do
!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving v_split_adv()'
   write(debug,*)
#endif
   return
   end subroutine v_split_adv

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
