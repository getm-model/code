#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  3D advection
!
! !INTERFACE:
   module advection_3d
!
! !DESCRIPTION:
!
!  This module do advection of scalars.  The module follows the same
!  convention as the other modules in 'getm'. The module is initialised
!  by calling 'init\_advection\_3d()'. In the time-loop 'do\_advection\_3d' is
!  called. 'do\_advection\_3d' is a wrapper routine which - dependent on the
!  actual advection scheme chosen - makes calls to the appropriate
!  subroutines, which may be done as one-step or multiple-step schemes.
!  The actual subroutines are coded in external FORTRAN files.
!  New advection schemes are easily implemented - at least from a program
!  point of view - since only this module needs to be changed.
!  Additional work arrays can easily be added following the stencil given
!  below. To add a new advection scheme three things must be done:
!
!  \begin{enumerate}
!  \item define
!  a unique constant to identify the scheme (see e.g.\ {\tt UPSTREAM}
!  and {\tt TVD})
!  \item adopt the {\tt select case} in {\tt do\_advection\_3d} and
!  \item  write the actual subroutine.
!  \end{enumerate}
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   use advection
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection_3d, do_advection_3d
!  Note (KK): hi and adv3d are used only in do_advection_3d
!             the other fields are provided for e.g. uv_advect_3d
#ifdef STATIC
   REALTYPE,dimension(I3DFIELD)               :: hi,adv3d
   REALTYPE,public,dimension(I3DFIELD)        :: uuadv,vvadv,huadv,hvadv,fadv3d
#else
   REALTYPE,dimension(:,:),allocatable        :: hi,adv3d
   REALTYPE,public,dimension(:,:),allocatable :: uuadv,vvadv,huadv,hvadv,fadv3d
#endif
   integer,public,parameter           :: HVSPLIT=3
   character(len=64),public,parameter :: adv_splits_3d(0:3) = &
      (/"no split: one 3D step",                              &
        "full step splitting: u + v + w",                     &
        "half step splitting: u/2 + v/2 + w + v/2 + u/2",     &
        "hor/ver splitting: uv + w"/)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_advection_3d
!
! !INTERFACE:
   subroutine init_advection_3d(method)
!
! !DESCRIPTION:
!
! Here, memory for some variables is allocated, which are then initialised to
! zero.
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection_3d() # ',Ncall
#endif

   LEVEL2 'init_advection_3d'

#ifndef STATIC
   allocate(hi(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (hi)'

   allocate(adv3d(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (adv3d)'

   allocate(uuadv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (uuadv)'

   allocate(vvadv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (vvadv)'

   allocate(huadv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (huadv)'

   allocate(hvadv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (hvadv)'

   allocate(fadv3d(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (fadv3d)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_advection_3d
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_advection_3d - 3D advection schemes \label{sec-do-advection-3d}
!
! !INTERFACE:
   subroutine do_advection_3d(dt,f,uu,vv,ww,hu,hv,ho,hn,                          &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                              dxu,dxv,dyu,dyv,arcd1,                              &
#endif
                              az,au,av,hor_adv,ver_adv,adv_split,AH,hires,advres)
!
! !DESCRIPTION:
!
! Here, advection terms for all three-dimensional state variables are
! calculated by means of a finite-volume approach (an exception
! is the possibility to directly calculate the momentum advection
! by a one-step three-dimensional upstream scheme,
! see {\tt uv\_advection\_3d}) and the advection step is carried out
! as a fractional advection time step. Those 3D variables may be defined on
! T-, U-, V- and W-points. The latter option is interesting for
! turbulent quantities, but is not coded yet. Inside this advection
! routine, it does not matter, wehre the advected variable is located
! on the grid. All finite volume fluxes and geometric coefficients
! need to be calculated before {\tt advection\_3d} is called.
!
! Originally, this 3D advection routine has been written for tracer
! equations. There,
! after multiplying the layer-integrated and transformed to
! curvilinear coordinates tracer equation (\ref{C_Layer_IntCurvi})
! with $mn$, the advective
! terms in this equation are discretised as follows.
!
! First advection term in (\ref{C_Layer_IntCurvi}):
! \begin{equation}\label{u_discr_advect}
! \left(mn\,\partial_{\cal X}\left(\frac{p_kc_k}{n}\right)\right)_{i,j}\approx
! \frac{
! p_{i,j,k}\tilde c^u_{i,j,k}\Delta y^u_{i,j}-
! p_{i-1,j,k}\tilde c^u_{i-1,j,k}\Delta y^u_{i-1,j}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}
! \end{equation}
!
! Second advection term in (\ref{C_Layer_IntCurvi}):
! \begin{equation}\label{v_discr_advect}
! \left(mn\,\partial_{\cal Y}\left(\frac{q_kc_k}{m}\right)\right)_{i,j}\approx
! \frac{
! q_{i,j,k}\tilde c^v_{i,j,k}\Delta y^v_{i,j}-
! q_{i,j-1,k}\tilde c^v_{i,j-1,k}\Delta y^v_{i,j-1}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}
! \end{equation}
!
! Vertical advective fluxes in (\ref{C_Layer_IntCurvi}):
! \begin{equation}\label{w_discr_advect}
! \left(\bar w_{k} \tilde c_{k}\right)_{i,j}\approx
! w_{i,j,k}\tilde c^w_{i,j,k}.
! \end{equation}
!
! The interfacial concentrations $\tilde c_{i,j,k}$ are calculated
! according to upwind or higher order directional split
! schemes, which are discussed in detail below and in sections
! \ref{sec-upstream-adv} - \ref{sec-fct-2dh-adv}.
!
! However, as said above, in the same way these routines may be applied
! to quantities on
! U-, V-, and W-points, if the transports and geometric coefficients
! are properly calculated.
!
! There are various combinations of advection schemes possible.
! The first selection is whether a one-step 3D first-order upstream method
! is cosen, or a fractional step method.
!
! The next selection is (if a fractional step method is selected)
! how to do the fractional steps (selection on {\tt adv\_split}). There
! are different options,
!
! \begin{enumerate}
! \item directional split with subsequent full steps in $x$-, $y$- and
! $z$-direction,
! \item split with subsequent half steps in $x$-, and $y$-direction, a
! full step in $z$-direction, and half steps in $y$- and $x$-direction.
! \item directional split into a 2D horizontal step and a 1D vertical step.
! \end{enumerate}
!
! For the 1D directional-split schemes, first-order upstream,
! ULTIMATE QUICKEST, and the Total Variation Diminishing (TVD) schemes
! Superbee, MUSCL, and P$_2$PDM are available.

! For the 2D horizontal step, an upstream scheme and a Flux-Corrected
! Transport (FCT) scheme have been coded.
!
! If the compiler option {\tt ITERATE\_VERT\_ADV} is chosen, the vertical
! advection is iterated as many times with reduced time step that
! the CFL criterium for vertical advection is fulfilled, see the routine
! {\tt w\_split\_it\_adv}.
!
! With the compiler option {\tt SLICE\_MODEL}, the advection in
! $y$-direction is not executed.
!
!
! !USES:
   use getm_timers, only: tic, toc, TIM_ADV3D, TIM_ADV3DH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                               :: dt,AH
   REALTYPE,dimension(I3DFIELD),intent(in)           :: uu,vv,ww,ho,hn,hu,hv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(I2DFIELD),intent(in)           :: dxu,dxv,dyu,dyv,arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)            :: az,au,av
   integer,intent(in)                                :: adv_split,hor_adv,ver_adv
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout)        :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(out),optional :: hires,advres
!
! !LOCAL VARIABLES:
   integer :: k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection_3d() # ',Ncall
#endif
   call tic(TIM_ADV3D)

   hi = ho
   adv3d = _ZERO_

    select case (adv_split)

      case(NOSPLIT)

         select case (hor_adv)

            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               do k=1,kmax
                  call adv_u_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   uu(:,:,k),ho(:,:,k),hu(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxu,dyu,arcd1,                      &
#endif
                                   au,_ONE_,hor_adv,az,AH,             &
                                   onestep_finalise=.false.)
               end do
#ifndef SLICE_MODEL
               do k=1,kmax
                  call adv_v_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   vv(:,:,k),ho(:,:,k),hv(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxv,dyv,arcd1,                      &
#endif
                                   av,_ONE_,hor_adv,az,AH,             &
                                   onestep_finalise=.false.)
               end do
#endif

            case (UPSTREAM_2DH)

               do k=1,kmax
                  call adv_upstream_2dh(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                        uu(:,:,k),vv(:,:,k),                &
                                        ho(:,:,k),hn(:,:,k),                &
                                        hu(:,:,k),hv(:,:,k),                &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                        dxv,dyu,dxu,dyv,arcd1,              &
#endif
                                        az,AH,onestep_finalise=.false.)
               end do

            case (FCT)

               do k=1,kmax
                  call adv_fct_2dh(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   uu(:,:,k),vv(:,:,k),                &
                                   ho(:,:,k),hn(:,:,k),                &
                                   hu(:,:,k),hv(:,:,k),                &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxv,dyu,dxu,dyv,arcd1,              &
#endif
                                   az,AH,onestep_finalise=.false.)
               end do

            case default

               FATAL 'hor_adv=',hor_adv,' is invalid'
               stop

         end select
         call adv_w_split_3d(dt,f,hi,adv3d,ww,az,_ONE_,ver_adv,itersmax, &
                             onestep_finalise=.true.)

      case(FULLSPLIT)

         do k=1,kmax
            call do_advection(dt,f(:,:,k),uu(:,:,k),vv(:,:,k),         &
                              hu(:,:,k),hv(:,:,k),ho(:,:,k),hn(:,:,k), &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                              dxu,dxv,dyu,dyv,arcd1,                   &
#endif
                              az,au,av,hor_adv,adv_split,AH,           &
                              Dires=hi(:,:,k),advres=adv3d(:,:,k))
         end do
         if (kmax .gt. 1) then
            call tic(TIM_ADV3DH)
            call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
            call wait_halo(D_TAG)
            call toc(TIM_ADV3DH)
            call adv_w_split_3d(dt,f,hi,adv3d,ww,az,_ONE_,ver_adv,itersmax)
         end if

      case(HALFSPLIT)

         select case (adv_scheme)

            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               do k=1,kmax
                  call adv_u_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   uu(:,:,k),ho(:,:,k),hu(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxu,dyu,arcd1,                      &
#endif
                                   au,_HALF_,hor_adv,az,AH)
               end do
#ifndef SLICE_MODEL
               call tic(TIM_ADV3DH)
               call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADV3DH)
               do k=1,kmax
                  call adv_v_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   vv(:,:,k),ho(:,:,k),hv(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxv,dyv,arcd1,                      &
#endif
                                   av,_HALF_,hor_adv,az,AH)
               end do
#endif
               if (kmax .gt. 1) then
                  call tic(TIM_ADV3DH)
                  call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
                  call wait_halo(D_TAG)
                  call toc(TIM_ADV3DH)
                  call adv_w_split_3d(dt,f,hi,adv3d,ww,az,_ONE_,ver_adv,itersmax)
               end if
#ifndef SLICE_MODEL
               call tic(TIM_ADV3DH)
               call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADV3DH)
               do k=1,kmax
                  call adv_v_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   vv(:,:,k),ho(:,:,k),hv(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxv,dyv,arcd1,                      &
#endif
                                   av,_HALF_,hor_adv,az,AH)
               end do
#endif
               call tic(TIM_ADV3DH)
               call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADV3DH)
               do k=1,kmax
                  call adv_u_split(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                   uu(:,:,k),ho(:,:,k),hu(:,:,k),      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                   dxu,dyu,arcd1,                      &
#endif
                                   au,_HALF_,hor_adv,az,AH)
               end do

            case((UPSTREAM_2DH),(FCT))

               stop 'do_advection_3d: hor_adv=',hor_adv,' not valid for adv_split=',adv_split

            case default

               stop 'do_advection_3d: hor_adv=',hor_adv,' is invalid'

         end select

      case(HVSPLIT)

         do k=1,kmax
            call do_advection(dt,f(:,:,k),uu(:,:,k),vv(:,:,k),         &
                              hu(:,:,k),hv(:,:,k),ho(:,:,k),hn(:,:,k), &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                              dxu,dxv,dyu,dyv,arcd1,                   &
#endif
                              az,au,av,hor_adv,NOSPLIT,AH,             &
                              Dires=hi(:,:,k),advres=adv3d(:,:,k))
         end do
         if (kmax .gt. 1) then
            call tic(TIM_ADV3DH)
            call update_3d_halo(f,f,az,imin,jmin,imax,jmax,kmax,D_TAG)
            call wait_halo(D_TAG)
            call toc(TIM_ADV3DH)
            call adv_w_split_3d(dt,f,hi,adv3d,ww,az,_ONE_,ver_adv,itersmax)
         end if

      case default

         stop 'do_advection_3d: adv_split=',adv_split,' is invalid'

   end select

   if (present(hires)) hires = hi
   if (present(advres)) advres = adv3d

   call toc(TIM_ADV3D)
#ifdef DEBUG
   write(debug,*) 'Leaving do_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_advection_3d
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  print_adv_settings_3d
!
! !INTERFACE:
   subroutine print_adv_settings_3d(adv_split,hor_adv,ver_adv,itersmax)
!
! !DESCRIPTION:
!
! Checks and prints out settings for 3D advection.
!
! !USES
   IMPLICIT NONE

! !INPUT PARAMETERS:
   integer,intent(in) :: adv_split,hor_adv,ver_adv,itersmax
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'print_adv_settings_3d() # ',Ncall
#endif

   select case (adv_split)
      case((NOSPLIT),(FULLSPLIT),(HALFSPLIT),(HVSPLIT))
      case default
         FATAL 'adv_split=',adv_split,' is invalid'
         stop
   end select

   select case (hor_adv)
      case((UPSTREAM),(UPSTREAM_2DH),(P2),(SUPERBEE),(MUSCL),(P2_PDM),(FCT))
      case default
         FATAL 'hor_adv=',hor_adv,' is invalid'
         stop
   end select

   select case (ver_adv)
      case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))
      case default
         FATAL 'ver_adv=',ver_adv,' is invalid'
         stop
   end select

   select case (adv_split)
      case((FULLSPLIT),(HALFSPLIT))
         select case (hor_adv)
            case((UPSTREAM_2DH),(FCT))
               FATAL 'hor_adv=',hor_adv,' not valid for adv_split=',adv_split
               stop
         end select
   end select

   LEVEL3 adv_splits_3d(adv_split)
   LEVEL3 ' horizontal: ',adv_schemes(hor_adv)
   LEVEL3 ' vertical  : ',adv_schemes(ver_adv)

   if (adv_split .eq. NOSPLIT) then
      LEVEL3 '             adv_split=',adv_split,' disables iteration
   else
      if (itersmax .gt. 1) then
         LEVEL3 '             with max ',itersmax,' iterations'
      else
         LEVEL3 '             without iteration'
      end if
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving print_adv_settings_3d()'
   write(debug,*)
#endif
   return
   end subroutine print_adv_settings_3d
!EOC
!-----------------------------------------------------------------------

   end module advection_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
