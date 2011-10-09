#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_advect - 2D advection of momentum \label{sec-uv-advect}
!
! !INTERFACE:
   subroutine uv_advect(U,V,DU,DV)
!
! !DESCRIPTION:
!
! The advective terms in the vertically integrated
! momentum equation are discretised in
! a momentum-conservative form. This is carried out here for the
! advective terms in the $U$-equation (\ref{UMOM}) and the
! $V$-equation (\ref{VMOM}) (after applying the curvilinear
! coordinate transformationand multiplying these
! equations with $mn$).
!
! First advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{U^2}{Dn}\right)\right)_{i,j}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(U_{i+1,j}+U_{i,j})\tilde u_{i+1,j}\Delta y^c_{i+1,j}-
! \frac12(U_{i,j}+U_{i-1,j})\tilde u_{i,j}\Delta y^c_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the inter-facial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i-1,j}}{D^u_{i-1,j}} & \mbox{ for } \frac12(U_{i,j}+U_{i-1,j})>0\\ \\
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}y\left(\frac{UV}{Dm}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(V_{i+1,j}+V_{i,j})\tilde u_{i,j}\Delta x^+_{i,j}-
! \frac12(V_{i+1,j-1}+V_{i,j-1})\tilde u_{i,j-1}\Delta x^+_{i,j-1}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the inter-facial
! velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ for } \frac12(V_{i+1,j,k}+V_{i,j,k})>0\\ \\
! \displaystyle
! \frac{U_{i,j+1}}{D^u_{i,j+1}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! First advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{UV}{Dn}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(U_{i,j+1}+U_{i,j})\tilde v_{i,j}\Delta y^+_{i,j}-
! \frac12(U_{i-1,j+1}+U_{i-1,j})\tilde v_{i-1,j}\Delta y^+_{i-1,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the interfacial
! velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ for } \frac12(U_{i+1,j}+U_{i,j})>0\\ \\
! \displaystyle
! \frac{V_{i+1,j}}{D^v_{i+1,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}\left(\frac{V^2}{Dm}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(V_{i,j+1}+V_{i,j})\tilde v_{i,j+1}\Delta x^c_{i,j+1}-
! \frac12(V_{i,j}+V_{i,j-1})\tilde v_{i,j}\Delta x^c_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the interfacial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j-1}}{D^v_{i,j-1}} & \mbox{ for } \frac12(V_{i,j}+V_{i,j-1})>0\\ \\
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxx,dyx,arud1,arvd1
#endif
   use m2d, only: dtm,vel_adv_split2d,vel_adv_scheme
   use variables_2d, only: UEx,VEx,fadv,Uadv,Vadv,DUadv,DVadv,maskadv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use variables_2d, only: dxadv,dyadv
#endif
   use advection, only: UPSTREAM,do_advection
   use halo_zones, only: update_2d_halo,wait_halo,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_UVADV,TIM_UVADVH
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: U,V,DU,DV
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect() # ',Ncall
#endif
   call tic(TIM_UVADV)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

!  Here begins dimensional split advection for u-velocity
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
!        the velocity to be transported
         fadv(i,j) = U(i,j)/DU(i,j)
         Uadv(i,j) = _HALF_*( U(i  ,j) + U(i+1,j) )
         Vadv(i,j) = _HALF_*( V(i  ,j) + V(i+1,j) )
!        Note (KK): DU only valid until imax+1
!                   therefore DUadv only valid until imax
         DUadv(i,j) = _HALF_*( DU(i  ,j) + DU(i+1,j) )
!        Note (KK): DV only valid until jmax+1
!                   therefore DVadv only valid until jmax+1
         DVadv(i,j) = _HALF_*( DV(i  ,j) + DV(i+1,j) )
         maskadv(i,j) = az(i+1,j)
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dxadv(i,j) = dxc(i+1,j)
         dyadv(i,j) = dyc(i+1,j)
#endif
      end do
   end do
!$OMP END DO
!$OMP MASTER
   if (vel_adv_scheme .ne. UPSTREAM) then
!     we need to update fadv(imax+HALO,jmin-HALO:jmax+HALO)
      call tic(TIM_UVADVH)
      call update_2d_halo(fadv,fadv,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)
      call toc(TIM_UVADVH)
   end if

   call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,DU,DU,        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                     dxadv,dxx,dyadv,dyx,arud1,                   &
#endif
                     au,maskadv,ax,                               &
                     vel_adv_scheme,vel_adv_split2d,_ZERO_,U_TAG, &
                     advres=UEx)
!$OMP END MASTER
!  OMP-NOTE: MASTER does not imply BARRIER
!$OMP BARRIER

!  Here begins dimensional split advection for v-velocity
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
!        the velocity to be transported
         fadv(i,j) = V(i,j)/DV(i,j)
         Uadv(i,j) = _HALF_*( U(i,j) + U(i,j+1) )
         Vadv(i,j) = _HALF_*( V(i,j) + V(i,j+1) )
!        Note (KK): DU only valid until imax+1
!                   therefore DUadv only valid until imax+1
         DUadv(i,j) = _HALF_*( DU(i,j) + DU(i,j+1) )
!        Note (KK): DV only valid until jmax+1
!                   therefore DVadv only valid until jmax
         DVadv(i,j) = _HALF_*( DV(i,j) + DV(i,j+1) )
         maskadv(i,j) = az(i,j+1)
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dxadv(i,j) = dxc(i,j+1)
         dyadv(i,j) = dyc(i,j+1)
#endif
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL

   if (vel_adv_scheme .ne. UPSTREAM) then
!     we need to update fadv(imin-HALO:imax+HALO,jmax+HALO)
      call tic(TIM_UVADVH)
      call update_2d_halo(fadv,fadv,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)
      call toc(TIM_UVADVH)
   end if

   call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,DV,DV,        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                     dxx,dxadv,dyx,dyadv,arvd1,                   &
#endif
                     av,ax,maskadv,                               &
                     vel_adv_scheme,vel_adv_split2d,_ZERO_,V_TAG, &
                     advres=VEx)

   call toc(TIM_UVADV)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
