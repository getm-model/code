#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: momentum - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine momentum(n,tausx,tausy,airp,ufirst)
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
! !INPUT/OUTPUT VARIABLES:
   logical, intent(inout)              :: ufirst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
!
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
   use domain, only: H,au,ax,min_depth,dry_u,Cori,coru
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,arud1,dyc,dxx
#else
   use domain, only: dx
#endif
   use variables_2d, only: UEuler
   use waves, only: waves_method,NO_WAVES
   use variables_waves, only: UStokes
   use domain, only: have_boundaries
   use variables_2d, only: dtm,D,z,UEx,U,DU,SlUx,Slru,ru,V,DV
   use bdy_2d, only: do_bdy_2d
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
   REALTYPE,dimension(E2DFIELD) :: work2d
   REALTYPE                  :: zp,zm,zx,tausu,Slr,Vloc,fV
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,zp,zm,zx,tausu,Slr,Vloc,fV,cord_curv)

!  calculate V at X-points (needed for Coriolis)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (ax(i,j) .ne. 0) then
#ifdef NEW_CORI
!           Espelid et al. [2000], IJNME 49, 1521-1545
            work2d(i,j) = _HALF_ * (  V(i  ,j)/sqrt(DV(i  ,j)) &
                                    + V(i+1,j)/sqrt(DV(i+1,j)) )
#else
            work2d(i,j) = _HALF_ * ( V(i,j) + V(i+1,j) )
#endif
         else
            work2d(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if ((au(i,j) .eq. 1) .or. (au(i,j) .eq. 2)) then
!           Semi-implicit treatment of Coriolis force
!           Note (KK): Between individual momentum routines no call
!                      to depth_update, therefore combination of new
!                      transport and old depth (impossible to fix).
!                      In the old implementation (fV calculated at the
!                      end of vmomentum), this was the default.
!                      In this new implementation at least the momentum
!                      called first gets consistent transports and depths.
            Vloc = _HALF_ * ( work2d(i,j-1) + work2d(i,j) )
#ifdef NEW_CORI
!           Espelid et al. [2000], IJNME 49, 1521-1545
            Vloc = Vloc * sqrt(DU(i,j))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(Vloc*(DYCIP1-DYC)-U(i,j)*(DXX-DXXJM1)) &
                      /DU(i,j)*ARUD1
            fV=(cord_curv+coru(i,j))*Vloc
#else
            fV=coru(i,j)*Vloc
#endif

            zp = max( z(i+1,j) , -H(i  ,j)+min( min_depth , D(i+1,j) ) )
            zm = max( z(i  ,j) , -H(i+1,j)+min( min_depth , D(i  ,j) ) )
            zx = ( zp - zm + (airp(i+1,j)-airp(i,j))*gammai ) / DXU
            tausu = _HALF_ * ( tausx(i,j) + tausx(i+1,j) )
            if (UEuler(i,j) .eq. _ZERO_) then
               Slr = Slru(i,j) / SMALL
            else
               Slr = Slru(i,j) / UEuler(i,j)
            end if
            Slr = max(_ZERO_,ru(i,j)/DU(i,j)+Slr)
            UEuler(i,j)=(UEuler(i,j)-dtm*(g*DU(i,j)*zx+dry_u(i,j)*&
                 (-tausu*rho_0i-fV+UEx(i,j)+SlUx(i,j))))/&
                 (_ONE_+dtm*Slr)
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL
! The rest of this sub is not easy to thread.

   if (have_boundaries) call do_bdy_2d(n,U_TAG)

#ifdef SLICE_MODEL
   j = jmax/2
   UEuler(imin:imax,j+1) = UEuler(imin:imax,j)
#endif

!  now u is calculated
   CALL tic(TIM_MOMENTUMH)
   call update_2d_halo(UEuler,UEuler,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)
   CALL toc(TIM_MOMENTUMH)
   call mirror_bdy_2d(UEuler,U_TAG)

   if (waves_method .ne. NO_WAVES) then
      U = UEuler + UStokes
   end if

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
   use domain, only: H,av,ax,min_depth,dry_v,Cori,corv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyv,arvd1,dyx,dxc
#else
   use domain, only: dy
#endif
   use variables_2d, only: VEuler
   use waves, only: waves_method,NO_WAVES
   use variables_waves, only: VStokes
   use domain, only: have_boundaries
   use variables_2d, only: dtm,D,z,VEx,V,DV,SlVx,Slrv,rv,U,DU
   use bdy_2d, only: do_bdy_2d
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
   REALTYPE,dimension(E2DFIELD) :: work2d
   REALTYPE                  :: zp,zm,zy,tausv,Slr,Uloc,fU
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,zp,zm,zy,tausv,Slr,Uloc,fU,cord_curv)

!  calculate U at X-points (needed for Coriolis)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (ax(i,j) .ne. 0) then
#ifdef NEW_CORI
!           Espelid et al. [2000], IJNME 49, 1521-1545
            work2d(i,j) = _HALF_ * (  U(i,j  )/sqrt(DU(i,j  )) &
                                    + U(i,j+1)/sqrt(DU(i,j+1)) )
#else
            work2d(i,j) = _HALF_ * ( U(i,j) + U(i,j+1) )
#endif
         else
            work2d(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if ((av(i,j) .eq. 1) .or. (av(i,j) .eq. 2)) then
!           Semi-implicit treatment of Coriolis force
!           Note (KK): Between individual momentum routines no call
!                      to depth_update, therefore combination of new
!                      transport and old depth (impossible to fix).
!                      In the old implementation (fV calculated at the
!                      end of vmomentum), this was the default.
!                      In this new implementation at least the momentum
!                      called first gets consistent transports and depths.
            Uloc = _HALF_ * ( work2d(i-1,j) + work2d(i,j) )
#ifdef NEW_CORI
!           Espelid et al. [2000], IJNME 49, 1521-1545
            Uloc = Uloc * sqrt(DV(i,j))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(V(i,j)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC)) &
                      /DV(i,j)*ARVD1
            fU=(cord_curv+corv(i,j))*Uloc
#else
            fU=corv(i,j)*Uloc
#endif

            zp = max( z(i,j+1) , -H(i,j  )+min( min_depth , D(i,j+1) ) )
            zm = max( z(i,j  ) , -H(i,j+1)+min( min_depth , D(i,j  ) ) )
            zy = ( zp - zm + (airp(i,j+1)-airp(i,j))*gammai ) / DYV
            tausv = _HALF_ * ( tausy(i,j) + tausy(i,j+1) )
            if (VEuler(i,j) .eq. _ZERO_) then
               Slr = Slrv(i,j) / SMALL
            else
               Slr = Slrv(i,j) / VEuler(i,j)
            end if
            Slr = max(_ZERO_,rv(i,j)/DV(i,j)+Slr)
            VEuler(i,j)=(VEuler(i,j)-dtm*(g*DV(i,j)*zy+dry_v(i,j)*&
                 (-tausv*rho_0i+fU+VEx(i,j)+SlVx(i,j))))/&
                 (_ONE_+dtm*Slr)
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL
! The rest of this sub is not easy to thread.

   if (have_boundaries) call do_bdy_2d(n,V_TAG)

#ifdef SLICE_MODEL
   j = jmax/2
   VEuler(imin:imax,j-1) = VEuler(imin:imax,j)
   VEuler(imin:imax,j+1) = VEuler(imin:imax,j)
#endif

!  now v is calculated
   CALL tic(TIM_MOMENTUMH)
   call update_2d_halo(VEuler,VEuler,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)
   CALL toc(TIM_MOMENTUMH)
   call mirror_bdy_2d(VEuler,V_TAG)

   if (waves_method .ne. NO_WAVES) then
      V = VEuler + VStokes
   end if

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
