#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_advect_3d - 3D momentum advection \label{sec-uv-advect-3d}
!
! !INTERFACE:
   subroutine uv_advect_3d()
!
! !DESCRIPTION:
!
! For the discretisation of the momentum advection terms, two
! conceptionally different methods have been implemented in GETM.
! The first is the straight-forward multidimensional advection
! scheme, which is here realised as the first-order upwind scheme,
! see paragraph {\bf Multidimensional approach} on page
! \pageref{uvadvect-multi}.
!
! In order to make use of the higher-order directional-split
! methods for tracers (see section \ref{sec-do-advection-3d}),
! an alternative method is implemented, in which the complete advection step
! is first made, and then the resulting advection terms,
! which are needed for the calculation of the slow terms, see equations
! (\ref{SxA}) and (\ref{SyA})) are calculated from this
! (see paragraph {\bf Directional-split approach} on page
! \pageref{uvadvect-direct}).
!
! The choice which of the two methods to be used is made by means
! of the compiler option {\tt UV\_TVD} which has to be set in the
! {\tt Makefile} of the application in order to activate the
! more accurate but computationally more demanding high-order directional-split
! method. The effect of hih-order advaction can be impressively studied
! by means of the freshwater lense test case described in detail by
! \cite{BURCHARDea02}.
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
!
! \paragraph{Multidimensional approach}\label{uvadvect-multi}
!
! The advective terms in the momentum equation are discretised in
! a momentum-conservative form. This is carried out here for the
! advective terms in the $u$-equation (\ref{uEqviCurvi}) and the
! $v$-equation (\ref{vEqviCurvi}) (after multiplying these
! equations with $mn$).
!
! First advection term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{u_kp_k}{n}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(p_{i+1,j,k}+p_{i,j,k})\tilde u_{i+1,j,k}\Delta y^c_{i+1,j}-
! \frac12(p_{i,j,k}+p_{i-1,j,k})\tilde u_{i,j,k}\Delta y^c_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For an upwind scheme, the inter-facial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i-1,j,k} & \mbox{ for } \frac12(p_{i,j,k}+p_{i-1,j,k})>0\\ \\
! u_{i,j,k} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}y\left(\frac{v_kp_k}{m}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(q_{i+1,j,k}+q_{i,j,k})\tilde u_{i,j,k}\Delta x^+_{i,j}-
! \frac12(q_{i+1,j-1,k}+q_{i,j-1,k})\tilde u_{i,j-1,k}\Delta x^+_{i,j-1}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For an upwind scheme, the inter-facial velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i,j,k} & \mbox{ for } \frac12(q_{i+1,j,k}+q_{i,j,k})>0\\ \\
! u_{i,j+1,k} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! First advection term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{v_kq_k}{n}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(p_{i,j+1,k}+p_{i,j,k})\tilde v_{i,j,k}\Delta y^+_{i,j}-
! \frac12(p_{i-1,j+1,k}+p_{i-1,j,k})\tilde v_{i-1,j,k}\Delta y^+_{i-1,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For an upwind scheme, the interfacial velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j,k} & \mbox{ for } \frac12(p_{i+1,j,k}+p_{i,j,k})>0\\ \\
! v_{i+1,j,k} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}\left(\frac{v_kq_k}{m}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(q_{i,j+1,k}+q_{i,j,k})\tilde v_{i,j+1,k}\Delta x^c_{i,j+1}-
! \frac12(q_{i,j,k}+q_{i,j-1,k})\tilde v_{i,j,k}\Delta x^c_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For an upwind scheme, the interfacial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j-1,k} & \mbox{ for } \frac12(q_{i,j,k}+q_{i,j-1,k})>0\\ \\
! v_{i,j,k} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! The vertical advection terms in equations (\ref{uEqviCurvi})
! and (\ref{vEqviCurvi})
! can be discretised in an upstream scheme as well.
!
! Vertical advective flux in (\ref{uEqviCurvi}):
! \begin{equation}
! \left(\bar w_{k} \tilde u_{k}\right)_{i,j}\approx
! \frac12(w_{i+1,j,k}+w_{i,j,k})\tilde u_{i,j,k}
! \end{equation}
!
! with
!
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i,j,k} & \mbox{ for } \frac12(w_{i+1,j,k}+w_{i,j,k}) > 0,\\ \\
! u_{i,j,k+1} & \mbox{ else.}
! \end{array}
! \right.
! \end{equation}
!
! Vertical advective flux in (\ref{vEqviCurvi}):
! \begin{equation}
! \left(\bar w_{k} \tilde v_{k}\right)_{i,j}\approx
! \frac12(w_{i,j+1,k}+w_{i,j,k})\tilde v_{i,j,k}
! \end{equation}
!
! with
!
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j,k} & \mbox{ for } \frac12(w_{i,j+1,k}+w_{i,j,k}) > 0,\\ \\
! v_{i,j,k+1} & \mbox{ else.}
! \end{array}
! \right.
! \end{equation}
!
!
! \paragraph{Directional-split approach}\label{uvadvect-direct}
!
! Multidimensional treatment of advective terms in three-dimensional
! models is often quite unhandy, especially when higher-order
! advection schemes are needed. On the other hand, directional-split
! methods (which update the advected fields
! in each directional step and then "forget" the advection terms)
! as discussed in section \ref{sec-do-advection-3d} on page
! \pageref{sec-do-advection-3d},
! cannot directly be used
! for momentum advection when the models are based on mode splitting as
! e.g.\ GETM. The reason for this is that the three-dimensional
! advection terms are also needed for calculating the slow terms of the
! barotropic (external) mode, see equations (\ref{SxA}) and (\ref{SyA}).
!
! The procedure suggested here is as follows. First, the pure momentum advection
! equations are formally solved with the directional-split method
! described in section \ref{sec-do-advection-3d}:
!
! \begin{equation}\label{uEqvi_adv}
! \begin{array}{l}
! \partial_t p_k
! +\partial_x(u_kp_k)+\partial_y(v_kp_k)
! +\bar w_k \tilde u_k -\bar w_{k-1} \tilde u_{k-1}
! = 0,
! \end{array}
! \end{equation}
!
! \begin{equation}\label{vEqvi_adv}
! \partial_t q_k
! +\partial_x(u_kq_k)+\partial_y(v_kq_k)
! +\bar w_k  \tilde v_k -\bar w_{k-1}  \tilde v_{k-1}
! =0.
! \end{equation}
!
! The new solutions $\hat p_{i,j,k}$ and $\hat q_{i,j,k}$
! are however not further used, but instead the resulting advective terms
! $-(\hat p_{i,j,k}-p_{i,j,k})/\Delta t$ and
! $-(\hat q_{i,j,k}-q_{i,j,k})/\Delta t$
! are later applied to the momentum equations (together with other
! processes such as horizontal diffusion, pressure gradients, etc.)
! and also used for the calculation of the slow terms in
! (\ref{SxA}) and (\ref{SyA}).
!
! With this method, all higher-order directional-split advection schemes
! are now available for the momentum advection. The advective
! fluxes needed for this have to be averaged from the conservative
! advective fluxes resulting from the continuity equation
! Continuity will
! still be retained due to the linearity of the continuity equation.
!
! !
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,au,av
   use m3d, only: vel_adv_split,vel_hor_adv,vel_ver_adv
   use variables_3d, only: dt,uu,vv,ww,hun,hvn,uuEx,vvEx
   use variables_3d, only: fadv3d,uuadv,vvadv,wwadv,huadv,hvadv
   use advection, only: UPSTREAM
   use advection_3d, only: do_advection_3d
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_UVADV3D,TIM_UVADV3DH
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect_3d() # ',Ncall
#endif
   call tic(TIM_UVADV3D)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

! Here begins dimensional split advection for u-velocity
   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
!           the velocity to be transported
            fadv3d(i,j,k) = uu(i,j,k)/hun(i,j,k)
            uuadv(i,j,k) = _HALF_*( uu(i,j,k) + uu(i+1,j,k) )
            vvadv(i,j,k) = _HALF_*( vv(i,j,k) + vv(i+1,j,k) )
            wwadv(i,j,k) = _HALF_*( ww(i,j,k) + ww(i+1,j,k) )
!           Note (KK): hun only valid until imax+1
!                      therefore huadv only valid until imax
            huadv(i,j,k) = _HALF_*( hun(i,j,k) + hun(i+1,j,k) )
!           Note (KK): hvn only valid until jmax+1
!                      therefore hvadv only valid until jmax+1
            hvadv(i,j,k) = _HALF_*( hvn(i,j,k) + hvn(i+1,j,k) )
         end do
      end do
!$OMP END DO NOWAIT
   end do

!$OMP BARRIER
!$OMP MASTER
   if (vel_hor_adv .ne. UPSTREAM) then
!     we need to update fadv3d(imax+HALO,jmin-HALO:jmax+HALO)
      call tic(TIM_UVADV3DH)
      call update_3d_halo(fadv3d,fadv3d,au,imin,jmin,imax,jmax,kmax,U_TAG)
      call wait_halo(U_TAG)
      call toc(TIM_UVADV3DH)
   end if

   call do_advection_3d(dt,fadv3d,uuadv,vvadv,wwadv,huadv,hvadv,hun,hun,    &
                        vel_hor_adv,vel_ver_adv,vel_adv_split,_ZERO_,U_TAG, &
                        advres=uuEx)
!$OMP END MASTER
!  OMP-NOTE: MASTER does not imply BARRIER
!$OMP BARRIER

!  Here begins dimensional split advection for v-velocity
   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
!           the velocity to be transported
            fadv3d(i,j,k) = vv(i,j,k)/hvn(i,j,k)
            uuadv(i,j,k) = _HALF_*( uu(i,j,k) + uu(i,j+1,k) )
            vvadv(i,j,k) = _HALF_*( vv(i,j,k) + vv(i,j+1,k) )
            wwadv(i,j,k) = _HALF_*( ww(i,j,k) + ww(i,j+1,k) )
!           Note (KK): hun only valid until imax+1
!                      therefore huadv only valid until imax+1
            huadv(i,j,k) = _HALF_*( hun(i,j,k) + hun(i,j+1,k) )
!           Note (KK): hvn only valid until jmax+1
!                      therefore hvadv only valid until jmax
            hvadv(i,j,k) = _HALF_*( hvn(i,j,k) + hvn(i,j+1,k) )
         end do
      end do
!$OMP END DO NOWAIT
   end do

!$OMP END PARALLEL

   if (vel_hor_adv .ne. UPSTREAM) then
!     we need to update fadv3d(imin-HALO:imax+HALO,jmax+HALO)
      call tic(TIM_UVADV3DH)
      call update_3d_halo(fadv3d,fadv3d,av,imin,jmin,imax,jmax,kmax,V_TAG)
      call wait_halo(V_TAG)
      call toc(TIM_UVADV3DH)
   end if

   call do_advection_3d(dt,fadv3d,uuadv,vvadv,wwadv,huadv,hvadv,hvn,hvn,    &
                        vel_hor_adv,vel_ver_adv,vel_adv_split,_ZERO_,V_TAG, &
                        advres=vvEx)

   call toc(TIM_UVADV3D)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
