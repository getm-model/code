#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: momentum - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine momentum(n,tausx,tausy,airp)
!
! !DESCRIPTION:
!
! This small routine calls the $U$-equation and the $V$-equation in an
! alternating sequence (UVVUUVVUUVVU), in order to provide higher
! accuracy and energy conservation for the explicit numerical treatment
! of the Coriolis term.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   ! For timer here: Only clock what is not taken at "next" level.
   use getm_timers, only: tic, toc, TIM_MOMENTUM
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
   REALTYPE, intent(in)                :: tausx(E2DFIELD)
   REALTYPE, intent(in)                :: tausy(E2DFIELD)
   REALTYPE, intent(in)                :: airp(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical                   :: ufirst=.false.
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'Momentum() # ',Ncall
#endif
   CALL tic(TIM_MOMENTUM)

   if(ufirst) then
      call umomentum(n,tausx,airp)
      call vmomentum(n,tausy,airp)
      ufirst = .false.
   else
      call vmomentum(n,tausy,airp)
      call umomentum(n,tausx,airp)
      ufirst = .true.
   end if

   CALL toc(TIM_MOMENTUM)
#ifdef DEBUG
   write(debug,*) 'Leaving momentum()'
   write(debug,*)
#endif
   return
   end subroutine momentum
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: umomentum - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine umomentum(n,tausx,airp)
!
! !DESCRIPTION:
!
! Here, the vertically integrated $U$-momentum equation (\ref{UMOM}) given
! on page \pageref{UMOM} including a
! number of slow terms is calculated. One slight modification is that
! for better stability of drying and flooding processes the slow friction
! term $S^x_F$ is now also multiplied with the parameter $\alpha$ defined
! in eq.\ (\ref{alpha}).
!
! Furthermore, the horizontal pressure gradient $\partial^*_x\zeta$ is
! modified in order to
! support drying and flooding, see figure \ref{figpressgrad} on page
! \pageref{figpressgrad} and the explanations in section \ref{Section_dry}.
! $\partial^*_x\zeta$ is now also considering the atmospheric pressure
! gradient at sea surface height.
!
! For numerical stability reasons, the $U$-momentum equation is here
! discretised in time such that the
! bed friction is treated explicitely:
!
!  \begin{equation}\label{Umom_discrete}
!  \displaystyle
!  U^{n+1}=\frac{U^n-\Delta t_m(gD\partial^*_x\zeta
!  +\alpha(-\frac{\tau_x^s}{\rho_0}-fV^n+U_{Ex}+S_A^x-S_D^x+S_B^x+S_F^x))}
!  {1+\Delta t_m\frac{R}{D^2}\sqrt{\left(U^n\right)^2+\left(V^n\right)^2}},
!  \end{equation}
!
!  with $U_{Ex}$ combining advection and diffusion of $U$, see routines
!  {\tt uv\_advect} (section \ref{sec-uv-advect} on page
!  \pageref{sec-uv-advect}) and {\tt uv\_diffusion}
!  (section \ref{sec-uv-diffusion} on page
!  \pageref{sec-uv-diffusion}). The slow terms
!  are calculated in the routine {\tt slow\_terms} documented in section
!  \ref{sec-slow-terms} on page \pageref{sec-slow-terms}.
!  In (\ref{Umom_discrete}), $U^{n+1}$ denotes the transport on the
!  new and $U^n$ and $V^n$ the transports on the old time level.
!
!  The Coriolis term $fU$ for the subsequent $V$-momentum is also calculated
!  here, by directly interpolating the $U$-transports to the V-points
!  or by a method suggested by \cite{ESPELIDea00} which takes the
!  varying water depths into account.
!
!  Some provisions for proper behaviour of the $U$-transports when
!  GETM runs as slice model are made as well, see section
!  \ref{Section_GETM_Slice} on page \pageref{Section_GETM_Slice}.
!
! !USES:
   use parameters, only: g,rho_0
   use domain, only: imin,imax,jmin,jmax
   use domain, only: H,au,av,min_depth,dry_u,Cori,corv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,arvd1,dxc,dyx
   use variables_2d, only: V
#else
   use domain, only: dx
#endif
   use variables_2d, only: dtm,D,z,UEx,U,DU,fV,SlUx,Slru,ru,fU,DV
   use bdy_2d, only: do_bdy_2d
   use m2d, only: have_boundaries
   use getm_timers,  only: tic, toc, TIM_MOMENTUMH
   use halo_zones, only : update_2d_halo,wait_halo,U_TAG
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
   REALTYPE, intent(in)                :: tausx(E2DFIELD),airp(E2DFIELD)
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer                   :: i,j
   REALTYPE                  :: zp,zm,zx,tausu,Slr,Uloc
   REALTYPE                  :: cord_curv=_ZERO_
   REALTYPE, save            :: gammai,rho_0i
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'umomentum() # ',Ncall
#endif

   if (first) then
      rho_0i = _ONE_ / rho_0
      gammai = _ONE_ / (rho_0*max(SMALL,g))
      first = .false.
   end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,zp,zm,zx,tausu,Slr,Uloc)

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if ((au(i,j) .eq. 1) .or. (au(i,j) .eq. 2)) then
            zp = max( z(i+1,j) , -H(i  ,j)+min( min_depth , D(i+1,j) ) )
            zm = max( z(i  ,j) , -H(i+1,j)+min( min_depth , D(i  ,j) ) )
            zx = ( zp - zm + (airp(i+1,j)-airp(i,j))*gammai ) / DXU
            tausu = _HALF_ * ( tausx(i,j) + tausx(i+1,j) )
            Slr = Slru(i,j)
!            if (U(i,j) .gt. _ZERO_) then
!               Slr = max( Slru(i,j) , _ZERO_ )
!            else
!               Slr = min( Slru(i,j) , _ZERO_ )
!            end if
            U(i,j)=(U(i,j)-dtm*(g*DU(i,j)*zx+dry_u(i,j)*&
                 (-tausu*rho_0i-fV(i,j)+UEx(i,j)+SlUx(i,j)+Slr)))/&
                 (_ONE_+dtm*ru(i,j)/DU(i,j))
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL
! The rest of this sub is not easy to thread.

   if (have_boundaries) call do_bdy_2d(n,U_TAG)

#ifdef SLICE_MODEL
   do i=imin,imax
      U(i,3)=U(i,2)
   end do
#endif

!  now u is calculated
   CALL tic(TIM_MOMENTUMH)
   call update_2d_halo(U,U,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)
   CALL toc(TIM_MOMENTUMH)
   call mirror_bdy_2d(U,U_TAG)

! Semi-implicit treatment of Coriolis force for V-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(av(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Uloc= &
             ( U(i,j  )/sqrt(DU(i,j  ))+ U(i-1,j  )/sqrt(DU(i-1,j  ))  &
             + U(i,j+1)/sqrt(DU(i,j+1))+ U(i-1,j+1)/sqrt(DU(i-1,j+1))) &
               *_QUART_*sqrt(DV(i,j))
#else
            Uloc=_QUART_*( U(i-1,j)+U(i,j)+U(i-1,j+1)+U(i,j+1))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(V(i,j)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC)) &
                      /DV(i,j)*ARVD1
            fU(i,j)=(cord_curv+corv(i,j))*Uloc
#else
            fU(i,j)=corv(i,j)*Uloc
#endif
         else
            fU(i,j)= _ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving umomentum()'
   write(debug,*)
#endif
   return
   end subroutine umomentum
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: vmomentum - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine vmomentum(n,tausy,airp)
!
! !DESCRIPTION:
!
! Here, the vertically integrated $V$-momentum equation (\ref{VMOM}) given
! on page \pageref{VMOM} including a
! number of slow terms is calculated. One slight modification is that
! for better stability of drying and flooding processes the slow friction
! term $S^y_F$ is now also multiplied with the parameter $\alpha$ defined
! in eq.\ (\ref{alpha}).
!
! Furthermore, the horizontal pressure gradient $\partial^*_y\zeta$ is
! modified in order to
! support drying and flooding, see figure \ref{figpressgrad} on page
! \pageref{figpressgrad} and the explanations in section \ref{Section_dry}.
! $\partial^*_y\zeta$ is now also considering the atmospheric pressure
! gradient at sea surface height.
!
! For numerical stability reasons, the $V$-momentum equation is here
! discretised in time such that the
! bed friction is treated explicitely:
!
!  \begin{equation}\label{Vmom_discrete}
!  \displaystyle
!  V^{n+1}=\frac{V^n-\Delta t_m(gD\partial^*_y\zeta
!  +\alpha(-\frac{\tau_y^s}{\rho_0}+fU^n+V_{Ex}+S_A^y-S_D^y+S_B^y+S_F^y))}
!  {1+\Delta t_m\frac{R}{D^2}\sqrt{\left(U^n\right)^2+\left(V^n\right)^2}},
!  \end{equation}
!
!  with $V_{Ex}$ combining advection and diffusion of $V$, see routines
!  {\tt uv\_advect} (section \ref{sec-uv-advect} on page
!  \pageref{sec-uv-advect}) and {\tt uv\_diffusion}
!  (section \ref{sec-uv-diffusion} on page
!  \pageref{sec-uv-diffusion}). The slow terms
!  are calculated in the routine {\tt slow\_terms} documented in section
!  \ref{sec-slow-terms} on page \pageref{sec-slow-terms}.
!  In (\ref{Vmom_discrete}), $V^{n+1}$ denotes the transport on the
!  new and $U^n$ and $V^n$ the transports on the old time level.
!
!  The Coriolis term $fV$ for the subsequent $U$-momentum is also calculated
!  here, by directly interpolating the $U$-transports to the U-points
!  or by a method suggested by \cite{ESPELIDea00} which takes the
!  varying water depths into account.
!
!  Some provisions for proper behaviour of the $V$-transports when
!  GETM runs as slice model are made as well, see section
!  \ref{Section_GETM_Slice} on page \pageref{Section_GETM_Slice}.
!
! !USES:
   use parameters, only: g,rho_0
   use domain, only: imin,imax,jmin,jmax
   use domain, only: H,au,av,min_depth,dry_v,Cori,coru
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyv,arud1,dxx,dyc
   use m2d, only: U
#else
   use domain, only: dy
#endif
   use variables_2d, only: dtm,D,z,VEx,V,DV,fU,SlVx,Slrv,rv,fV,DU
   use bdy_2d, only: do_bdy_2d
   use m2d, only: have_boundaries
   use getm_timers,  only: tic, toc, TIM_MOMENTUMH
   use halo_zones, only : update_2d_halo,wait_halo,V_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
   REALTYPE, intent(in)                :: tausy(E2DFIELD),airp(E2DFIELD)
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer                   :: i,j
   REALTYPE                  :: zp,zm,zy,tausv,Slr,Vloc
   REALTYPE                  :: cord_curv=_ZERO_
   REALTYPE, save            :: gammai,rho_0i
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'vmomentum() # ',Ncall
#endif

   if (first) then
      rho_0i = _ONE_ / rho_0
      gammai = _ONE_ / (rho_0*max(SMALL,g))
      first = .false.
   end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,zp,zm,zy,tausv,Slr,Vloc)

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if ((av(i,j) .eq. 1) .or. (av(i,j) .eq. 2)) then
            zp = max( z(i,j+1) , -H(i,j  )+min( min_depth , D(i,j+1) ) )
            zm = max( z(i,j  ) , -H(i,j+1)+min( min_depth , D(i,j  ) ) )
            zy = ( zp - zm + (airp(i,j+1)-airp(i,j))*gammai ) / DYV
            tausv = _HALF_ * ( tausy(i,j) + tausy(i,j+1) )
            Slr = Slrv(i,j)
!            if (V(i,j) .gt. _ZERO_) then
!               Slr = max( Slrv(i,j) , _ZERO_ )
!            else
!               Slr = min( Slrv(i,j) , _ZERO_ )
!            end if
            V(i,j)=(V(i,j)-dtm*(g*DV(i,j)*zy+dry_v(i,j)*&
                 (-tausv*rho_0i+fU(i,j)+VEx(i,j)+SlVx(i,j)+Slr)))/&
                 (_ONE_+dtm*rv(i,j)/DV(i,j))
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL
! The rest of this sub is not easy to thread.

   if (have_boundaries) call do_bdy_2d(n,V_TAG)

#ifdef SLICE_MODEL
   do i=imin,imax
      V(i,1)=V(i,2)
      V(i,3)=V(i,2)
   end do
#endif

!  now v is calculated
   CALL tic(TIM_MOMENTUMH)
   call update_2d_halo(V,V,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)
   CALL toc(TIM_MOMENTUMH)
   call mirror_bdy_2d(V,V_TAG)

!  Semi-implicit treatment of Coriolis force for U-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(au(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Vloc = &
            ( V(i,j  )/sqrt(DV(i,j  ))+ V(i+1,j  )/sqrt(DV(i+1,j  )) + &
              V(i,j-1)/sqrt(DV(i,j-1))+ V(i+1,j-1)/sqrt(DV(i+1,j-1)))  &
              *_QUART_*sqrt(DU(i,j))
#else
            Vloc = _QUART_*( V(i,j-1)+ V(i+1,j-1)+V(i,j)+V(i+1,j))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(Vloc*(DYCIP1-DYC)-U(i,j)*(DXX-DXXJM1)) &
                       /DU(i,j)*ARUD1
            fV(i,j)=(cord_curv+coru(i,j))*Vloc
#else
            fV(i,j)=coru(i,j)*Vloc
#endif
         else
            fV(i,j) = _ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving vmomentum()'
   write(debug,*)
#endif
   return
   end subroutine vmomentum
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
