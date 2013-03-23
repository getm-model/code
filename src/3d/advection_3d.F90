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
!  This module does 3D advection of scalars. The module follows the same
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
   use advection
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection_3d, do_advection_3d,print_adv_settings_3d
   integer,public                        :: adv_ver_iterations=1
   integer,public,parameter              :: HVSPLIT=3,W_TAG=33
   character(len=64),public,parameter    :: adv_splits_3d(0:3) = &
            (/"no split: one 3D uvw step                     ",  &
              "full step splitting: u + v + w                ",  &
              "half step splitting: u/2 + v/2 + w + v/2 + u/2",  &
              "hor/ver splitting: uv + w                     "/)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine adv_split_w(dt,f,hi,adv3d,ww,ho,              &
                             splitfac,scheme,action,tag_3d,az, &
                             mask_finalise)
         use domain, only: imin,imax,jmin,jmax,kmax
         IMPLICIT NONE
         REALTYPE,intent(in)                             :: dt,splitfac
         REALTYPE,dimension(I3DFIELD),intent(in)         :: ww,ho
         integer,intent(in)                              :: scheme,action,tag_3d
         integer,dimension(E2DFIELD),intent(in)          :: az
         logical,dimension(E2DFIELD),intent(in),optional :: mask_finalise
         REALTYPE,dimension(I3DFIELD),intent(inout)      :: f,hi,adv3d
      end subroutine adv_split_w
   end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_advection_3d
!
! !INTERFACE:
   subroutine init_advection_3d()
!
! !DESCRIPTION:
!
! Allocates memory.
!
! !USES:
   IMPLICIT NONE
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection_3d() # ',Ncall
#endif

   LEVEL2 'init_advection_3d'

#ifdef ITERATE_VERT_ADV
   if (adv_ver_iterations .eq. 1) then
      adv_ver_iterations = 200
      LEVEL3 'changed number of maximum iterations'
      LEVEL3 'from 1 to 200 because of obsolete'
      LEVEL3 'compiler option -DITERATE_VERT_ADV'
   end if
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
   subroutine do_advection_3d(dt,f,uu,vv,ww,hu,hv,ho,hn,       &
                              split,hscheme,vscheme,AH,tag_3d, &
                              hires,advres)
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
! turbulent quantities. Inside this advection
! routine, it does not matter, wehre the advected variable is located
! on the grid. All finite volume fluxes and geometric coefficients
! need to be calculated before {\tt do\_advection\_3d} is called.
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
! \ref{sec-do-advection} and \ref{sec-w-split-adv}.
!
! However, as said above, in the same way these routines may be applied
! to quantities on
! U-, V-, and W-points, if the transports are properly calculated.
!
! There are various combinations of advection schemes possible.
!
! The options for {\tt split} are:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt split = NOSPLIT}: & no split (one 3D uvw step) \\
! {\tt split = FULLSPLIT}: & full step splitting (u + v + w) \\
! {\tt split = HALFSPLIT}: & half step splitting (u/2 + v/2 + w + v/2 + u/2) \\
! {\tt split = HVSPLIT}: & hor./ver. splitting (uv + w) \\
! \end{tabular}
!
! \vspace{0.5cm}
!
! The options for the horizontal scheme {\tt hscheme} are:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt scheme = NOADV}: & advection disabled \\
! {\tt scheme = UPSTREAM}: & first-order upstream (monotone) \\
! {\tt scheme = UPSTREAM\_2DH}: & 2DH upstream with forced monotonicity \\
! {\tt scheme = P2}: & third-order polynomial (non-monotone) \\
! {\tt scheme = SUPERBEE}: & second-order TVD (monotone) \\
! {\tt scheme = MUSCL}: & second-order TVD (monotone) \\
! {\tt scheme = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! {\tt scheme = J7}: & 2DH Arakawa J7 \\
! {\tt scheme = FCT}: & 2DH FCT with forced monotonicity \\
! {\tt scheme = P2\_2DH}: & 2DH P2 with forced monotonicity \\
! \end{tabular}
!
! \vspace{0.5cm}
!
! The options for the vertical scheme {\tt vscheme} are:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt scheme = NOADV}: & advection disabled \\
! {\tt scheme = UPSTREAM}: & first-order upstream (monotone) \\
! {\tt scheme = P2}: & third-order polynomial (non-monotone) \\
! {\tt scheme = SUPERBEE}: & second-order TVD (monotone) \\
! {\tt scheme = MUSCL}: & second-order TVD (monotone) \\
! {\tt scheme = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! \end{tabular}
!
! \vspace{0.5cm}
!
! With the compiler option {\tt SLICE\_MODEL}, the advection in
! meridional direction is not executed.
!
!
! !USES:
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_ADV3D,TIM_ADV3DH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                               :: dt,AH
   REALTYPE,dimension(I3DFIELD),intent(in)           :: uu,vv,ww,ho,hn,hu,hv
   integer,intent(in)                                :: split,hscheme,vscheme,tag_3d
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout)        :: f
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(out),optional :: hires,advres
!
! !LOCAL VARIABLES:
   type(t_adv_grid),pointer     :: adv_grid
   REALTYPE,dimension(I3DFIELD) :: hi,adv3d
   integer                      :: tag,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection_3d() # ',Ncall
#endif
   call tic(TIM_ADV3D)

   select case (tag_3d)
      case(H_TAG,D_TAG)
         tag = H_TAG
         adv_grid => adv_gridH
      case(U_TAG)
         tag = U_TAG
         adv_grid => adv_gridU
      case(V_TAG)
         tag = V_TAG
         adv_grid => adv_gridV
      case(W_TAG)
         tag = H_TAG
         adv_grid => adv_gridH
      case default
         stop 'do_advection_3d: tag_3d is invalid'
   end select

   hi = ho
   adv3d = _ZERO_

   if (hscheme.ne.NOADV .or. vscheme.ne.NOADV) then

      if (hscheme .eq. NOADV) then
         if (kmax .gt. 1) then
            call adv_split_w(dt,f,hi,adv3d,ww,ho,                           &
                             _ONE_,vscheme,SPLIT_UPDATE,tag_3d,adv_grid%az)
         end if
      else if (vscheme .eq. NOADV) then
         do k=1,kmax
            call do_advection(dt,f(:,:,k),uu(:,:,k),vv(:,:,k),         &
                              hu(:,:,k),hv(:,:,k),ho(:,:,k),hn(:,:,k), &
                              split,hscheme,AH,tag,                    &
                              Dires=hi(:,:,k),advres=adv3d(:,:,k))
         end do
      else

         select case (split)

            case(NOSPLIT)

               select case (hscheme)

                  case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

                     do k=1,kmax
                        call adv_split_u(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         uu(:,:,k),ho(:,:,k),hu(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                         _ONE_,hscheme,NOSPLIT_NOFINALISE,AH,       &
                                         adv_grid%mask_uflux,adv_grid%mask_uupdate)
#ifndef SLICE_MODEL
                        call adv_split_v(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         vv(:,:,k),ho(:,:,k),hv(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                         _ONE_,hscheme,NOSPLIT_NOFINALISE,AH,       &
                                         adv_grid%mask_vflux,adv_grid%mask_vupdate)
#endif
                     end do

                  case (UPSTREAM_2DH)

                     do k=1,kmax
                        call adv_upstream_2dh(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                              uu(:,:,k),vv(:,:,k),                &
                                              ho(:,:,k),hn(:,:,k),                &
                                              hu(:,:,k),hv(:,:,k),                &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                              adv_grid%dxv,adv_grid%dyu,          &
                                              adv_grid%dxu,adv_grid%dyv,          &
                                              adv_grid%arcd1,                     &
#endif
                                              NOSPLIT_NOFINALISE,AH,adv_grid%az)
                     end do

                  case(J7)

                     do k=1,kmax
                        call adv_arakawa_j7_2dh(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                                uu(:,:,k),vv(:,:,k),                &
                                                ho(:,:,k),hn(:,:,k),                &
                                                hu(:,:,k),hv(:,:,k),                &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                                adv_grid%dxv,adv_grid%dyu,          &
                                                adv_grid%dxu,adv_grid%dyv,          &
                                                adv_grid%arcd1,                     &
#endif
                                                NOSPLIT_NOFINALISE,AH,adv_grid%az,  &
                                                adv_grid%mask_uflux,                &
                                                adv_grid%mask_vflux,                &
                                                adv_grid%mask_xflux)
                     end do

                  case (FCT)

                     do k=1,kmax
                        call adv_fct_2dh(.true.,dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                         uu(:,:,k),vv(:,:,k),                       &
                                         ho(:,:,k),hn(:,:,k),                       &
                                         hu(:,:,k),hv(:,:,k),                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxv,adv_grid%dyu,                 &
                                         adv_grid%dxu,adv_grid%dyv,                 &
                                         adv_grid%arcd1,                            &
#endif
                                         NOSPLIT_NOFINALISE,AH,adv_grid%az,         &
                                         adv_grid%mask_uflux,                       &
                                         adv_grid%mask_vflux)
                     end do

                  case (P2_2DH)

                     do k=1,kmax
                        call adv_fct_2dh(.false.,dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k), &
                                         uu(:,:,k),vv(:,:,k),                        &
                                         ho(:,:,k),hn(:,:,k),                        &
                                         hu(:,:,k),hv(:,:,k),                        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxv,adv_grid%dyu,                  &
                                         adv_grid%dxu,adv_grid%dyv,                  &
                                         adv_grid%arcd1,                             &
#endif
                                         NOSPLIT_NOFINALISE,AH,adv_grid%az,          &
                                         adv_grid%mask_uflux,                        &
                                         adv_grid%mask_vflux)
                     end do

                  case default

                     stop 'do_advection_3d: hscheme is invalid'

               end select
!              Note (KK): here adv_split_w must be called in any case !!!
               call adv_split_w(dt,f,hi,adv3d,ww,ho,                   &
                                _ONE_,vscheme,NOSPLIT_FINALISE,tag_3d, &
                                adv_grid%az,                           &
                                mask_finalise=adv_grid%mask_finalise)

            case(FULLSPLIT)

               do k=1,kmax
                  call do_advection(dt,f(:,:,k),uu(:,:,k),vv(:,:,k),         &
                                    hu(:,:,k),hv(:,:,k),ho(:,:,k),hn(:,:,k), &
                                    FULLSPLIT,hscheme,AH,tag,                &
                                    Dires=hi(:,:,k),advres=adv3d(:,:,k))
               end do
               if (kmax .gt. 1) then
                  call adv_split_w(dt,f,hi,adv3d,ww,ho,               &
                                   _ONE_,vscheme,SPLIT_UPDATE,tag_3d, &
                                   adv_grid%az)
               end if

            case(HALFSPLIT)

               select case (hscheme)

                  case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

                     do k=1,kmax
                        call adv_split_u(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         uu(:,:,k),ho(:,:,k),hu(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                         _HALF_,hscheme,SPLIT_UPDATE,AH,            &
                                         adv_grid%mask_uflux,adv_grid%mask_uupdate)
                     end do
#ifndef SLICE_MODEL
#ifdef GETM_PARALLEL
                     if (hscheme.ne.UPSTREAM .and. tag_3d.eq.V_TAG) then
!                       we need to update f(imin:imax,jmax+HALO)
                        call tic(TIM_ADV3DH)
                        call update_3d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,kmax,D_TAG)
                        call wait_halo(D_TAG)
                        call toc(TIM_ADV3DH)
                     end if
#endif
                     do k=1,kmax
                        call adv_split_v(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         vv(:,:,k),ho(:,:,k),hv(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                         _HALF_,hscheme,SPLIT_UPDATE,AH,            &
                                         adv_grid%mask_vflux,adv_grid%mask_vupdate)
                     end do
#endif

                     if (kmax .gt. 1) then
                        call adv_split_w(dt,f,hi,adv3d,ww,ho,               &
                                         _ONE_,vscheme,SPLIT_UPDATE,tag_3d, &
                                         adv_grid%az)
                     end if

#ifndef SLICE_MODEL
#ifdef GETM_PARALLEL
                     if (hscheme .eq. UPSTREAM) then
                        if (tag_3d .eq. V_TAG) then
!                          we need to update f(imin-1:imax+1,jmax+1)
!                          KK-TODO: if external hv was halo-updated this halo-update is not necessary
                           call tic(TIM_ADV3DH)
                           call update_3d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,kmax,D_TAG)
                           call wait_halo(D_TAG)
                           call toc(TIM_ADV3DH)
                        end if
                     else
!                       we need to update f(imin:imax,jmin-HALO:jmin-1)
!                       we need to update f(imin:imax,jmax+1:jmax+HALO)
                        call tic(TIM_ADV3DH)
                        call update_3d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,kmax,D_TAG)
                        call wait_halo(D_TAG)
                        call toc(TIM_ADV3DH)
                     end if
#endif
                     do k=1,kmax
                        call adv_split_v(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         vv(:,:,k),ho(:,:,k),hv(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                         _HALF_,hscheme,SPLIT_UPDATE,AH,            &
                                         adv_grid%mask_vflux,adv_grid%mask_vupdate)
                     end do
#endif
#ifdef GETM_PARALLEL
                     if (hscheme .eq. UPSTREAM) then
                        if (tag_3d .eq. U_TAG) then
!                          we need to update f(imax+1,jmin:jmax)
!                          KK-TODO: if external hu was halo-updated this halo-update is not necessary
                           call tic(TIM_ADV3DH)
                           call update_3d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,kmax,D_TAG)
                           call wait_halo(D_TAG)
                           call toc(TIM_ADV3DH)
                        end if
                     else
!                       we need to update f(imin-HALO:imin-1,jmin:jmax)
!                       we need to update f(imax+1:imax+HALO,jmin:jmax)
                        call tic(TIM_ADV3DH)
                        call update_3d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,kmax,D_TAG)
                        call wait_halo(D_TAG)
                        call toc(TIM_ADV3DH)
                     end if
#endif
                     do k=1,kmax
                        call adv_split_u(dt,f(:,:,k),hi(:,:,k),adv3d(:,:,k),        &
                                         uu(:,:,k),ho(:,:,k),hu(:,:,k),             &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                         adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                         _HALF_,hscheme,SPLIT_UPDATE,AH,            &
                                         adv_grid%mask_uflux,adv_grid%mask_uupdate)
                     end do

                  case((UPSTREAM_2DH),(J7),(FCT),(P2_2DH))

                     stop 'do_advection_3d: hscheme not valid for split'

                  case default

                     stop 'do_advection_3d: hscheme is invalid'

               end select

            case(HVSPLIT)

               do k=1,kmax
                  call do_advection(dt,f(:,:,k),uu(:,:,k),vv(:,:,k),         &
                                    hu(:,:,k),hv(:,:,k),ho(:,:,k),hn(:,:,k), &
                                    NOSPLIT,hscheme,AH,tag,                  &
                                    Dires=hi(:,:,k),advres=adv3d(:,:,k))
               end do
               if (kmax .gt. 1) then
                  call adv_split_w(dt,f,hi,adv3d,ww,ho,               &
                                   _ONE_,vscheme,SPLIT_UPDATE,tag_3d, &
                                   adv_grid%az)
               end if

            case default

               stop 'do_advection_3d: split is invalid'

         end select

      end if

#ifdef SLICE_MODEL
      j = jmax/2
      f(:,j+1,:)     = f(:,j,:)
      hi(:,j+1,:)    = hi(:,j,:)
      adv3d(:,j+1,:) = adv3d(:,j,:)
      if (tag_3d .eq. V_TAG) then
         f(:,j-1,:)     = f(:,j,:)
         hi(:,j-1,:)    = hi(:,j,:)
         adv3d(:,j-1,:) = adv3d(:,j,:)
      end if
#endif

   end if

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
   subroutine print_adv_settings_3d(split,hscheme,vscheme,AH)
!
! !DESCRIPTION:
!
! Checks and prints out settings for 3D advection.
!
! !USES
   IMPLICIT NONE

! !INPUT PARAMETERS:
   integer,intent(inout):: split
   integer,intent(in)   :: hscheme,vscheme
   REALTYPE,intent(in)  :: AH
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

   if (hscheme.ne.NOADV .or. vscheme.ne.NOADV) then
      select case (split)
         case((NOSPLIT),(FULLSPLIT),(HALFSPLIT),(HVSPLIT))
         case default
            FATAL 'adv_split=',split,' is invalid'
            stop
      end select
   end if

   select case (hscheme)
      case((NOADV),(UPSTREAM),(UPSTREAM_2DH),(P2),(SUPERBEE),(MUSCL),(P2_PDM),(J7),(FCT),(P2_2DH))
      case default
         FATAL 'hor_adv=',hscheme,' is invalid'
         stop
   end select

   select case (vscheme)
      case((NOADV),(UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))
      case default
         FATAL 'ver_adv=',vscheme,' is invalid'
         stop
   end select

   if (vscheme .eq. NOADV) then
      if (split .eq. HVSPLIT) split = NOSPLIT
   end if

   if (hscheme .ne. NOADV) then
      select case (split)
         case((FULLSPLIT),(HALFSPLIT))
            select case (hscheme)
               case((UPSTREAM_2DH),(J7),(FCT),(P2_2DH))
                  FATAL 'hor_adv=',hscheme,' not valid for adv_split=',split
                  stop
            end select
      end select
      if (vscheme .eq. NOADV) then
         LEVEL3 trim(adv_splits(split))
      else
         LEVEL3 trim(adv_splits_3d(split))
      end if
   end if

   LEVEL3 ' horizontal: ',trim(adv_schemes(hscheme))

   if (hscheme .ne. NOADV) then
      if (AH .gt. _ZERO_) then
         LEVEL3 '             with AH=',AH
      else
         LEVEL3 '             without diffusion'
      end if
   end if

   LEVEL3 ' vertical  : ',trim(adv_schemes(vscheme))

   if (vscheme .ne. NOADV) then
      if (split .eq. NOSPLIT) then
         LEVEL3 '             adv_split=',split,' disables iteration'
      else
         if (adv_ver_iterations .gt. 1) then
            LEVEL3 '             with max ',adv_ver_iterations,' iterations'
         else
            LEVEL3 '             without iteration'
         end if
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
