#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  2D advection
!
! !INTERFACE:
   module advection
!
! !DESCRIPTION:
!
!  This module do advection of scalars.  The module follows the same
!  convention as the other modules in 'getm'. The module is initialised
!  by calling 'init\_advection()'. In the time-loop 'do\_advection()' is
!  called. 'do\_advection' is a wrapper routine which - dependent on the
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
!  \item adopt the {\tt select case} in {\tt do\_advection} and
!  \item  write the actual subroutine.
!  \end{enumerate}
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   IMPLICIT NONE

   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection,do_advection,print_adv_settings
   public adv_u_split,adv_v_split,adv_upstream_2dh,adv_fct_2dh

   type, public :: t_adv_grid
      logical,dimension(:,:),pointer :: mask_uflux,mask_uupdate
      logical,dimension(:,:),pointer :: mask_vflux,mask_vupdate
      logical,dimension(:,:),pointer :: mask_finalise
      integer,dimension(:,:),pointer :: az
#if defined(SPHERICAL) || defined(CURVILINEAR)
      REALTYPE,dimension(:,:),pointer :: dxu,dyu,dxv,dyv,arcd1
#endif
   end type t_adv_grid

   type(t_adv_grid),public,target :: adv_gridH,adv_gridU,adv_gridV

#ifdef STATIC
   logical,dimension(E2DFIELD),target         :: mask_updateH
   logical,dimension(E2DFIELD),target         :: mask_uupdateU,mask_uupdateV
   logical,dimension(E2DFIELD),target         :: mask_vupdateU,mask_vupdateV
   logical,dimension(_IRANGE_HALO_-1,_JRANGE_HALO_),target :: mask_ufluxV
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),target :: mask_vfluxU
   REALTYPE,public,dimension(E2DFIELD)        :: flux
   REALTYPE,dimension(E2DFIELD)               :: Di,adv
#else
   logical,dimension(:,:),allocatable,target  :: mask_updateH
   logical,dimension(:,:),allocatable,target  :: mask_uupdateU,mask_uupdateV
   logical,dimension(:,:),allocatable,target  :: mask_vupdateU,mask_vupdateV
   logical,dimension(:,:),allocatable,target  :: mask_ufluxV,mask_vfluxU
   REALTYPE,public,dimension(:,:),allocatable :: flux
   REALTYPE,dimension(:,:),allocatable        :: Di,adv
#endif
   integer,public,parameter           :: NOSPLIT=0,FULLSPLIT=1,HALFSPLIT=2
   character(len=64),public,parameter :: adv_splits(0:2) = &
                  (/"no split: one 2D uv step",            &
                    "full step splitting: u + v",          &
                    "half step splitting: u/2 + v + u/2"/)
   integer,public,parameter           :: UPSTREAM=1,UPSTREAM_2DH=2,P2=3
   integer,public,parameter           :: SUPERBEE=4,MUSCL=5,P2_PDM=6,FCT=7
   character(len=64),public,parameter :: adv_schemes(7) =  &
      (/"upstream advection (first-order, monotone)",      &
        "2DH-upstream advection with forced monotonicity", &
        "P2 advection (third-order, non-monotone)",        &
        "TVD-Superbee advection (second-order, monotone)", &
        "TVD-MUSCL advection (second-order, monotone)",    &
        "TVD-P2-PDM advection (third-order, monotone)",    &
        "2DH-FCT advection"/)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine adv_u_split(dt,f,Di,adv,U,Do,DU,            &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxu,dyu,arcd1,                  &
#endif
                             splitfac,scheme,AH,             &
                             mask_flux,mask_update,          &
                             nosplit_finalise,mask_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                             :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)         :: U,Do,DU
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(_IRANGE_HALO_-1,_JRANGE_HALO_),intent(in) :: dxu,dyu
         REALTYPE,dimension(E2DFIELD),intent(in)         :: arcd1
#endif
         integer,intent(in)                              :: scheme
         logical,dimension(_IRANGE_HALO_-1,_JRANGE_HALO_),intent(in) :: mask_flux
         logical,dimension(E2DFIELD),intent(in)          :: mask_update
         logical,intent(in),optional                     :: nosplit_finalise
         logical,dimension(E2DFIELD),intent(in),optional :: mask_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout)      :: f,Di,adv
      end subroutine adv_u_split

      subroutine adv_v_split(dt,f,Di,adv,V,Do,DV,            &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyv,arcd1,                  &
#endif
                             splitfac,scheme,AH,             &
                             mask_flux,mask_update,          &
                             nosplit_finalise,mask_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                             :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)         :: V,Do,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
         REALTYPE,dimension(E2DFIELD),intent(in)         :: arcd1
#endif
         integer,intent(in)                              :: scheme
         logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: mask_flux
         logical,dimension(E2DFIELD),intent(in)          :: mask_update
         logical,intent(in),optional                     :: nosplit_finalise
         logical,dimension(E2DFIELD),intent(in),optional :: mask_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout)      :: f,Di,adv
      end subroutine adv_v_split

      subroutine adv_upstream_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                  dxv,dyu,dxu,dyv,arcd1,       &
#endif
                                  az,AH,nosplit_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(_IRANGE_HALO_-1,_JRANGE_HALO_),intent(in) :: dxu,dyu
         REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
         REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         logical,intent(in),optional                :: nosplit_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_upstream_2dh

      subroutine adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyu,dxu,dyv,arcd1,       &
#endif
                             az,AH,nosplit_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(_IRANGE_HALO_-1,_JRANGE_HALO_),intent(in) :: dxu,dyu
         REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
         REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         logical,intent(in),optional                :: nosplit_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_fct_2dh
   end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_advection
!
! !INTERFACE:
   subroutine init_advection()
!
! !DESCRIvPTION:
!
! Here, memory for some variables is allocated, which are then initialised to
! zero.
!
! !USES
   use domain, only: az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,arcd1,dxu,dyu,arud1,dxv,dyv,arvd1,dxx,dyx
#endif
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection() # ',Ncall
#endif

   LEVEL2 'init_advection'

#ifndef STATIC
   allocate(mask_updateH(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_updateH)'

   allocate(mask_uupdateU(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_uupdateU)'

   allocate(mask_uupdateV(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_uupdateV)'

   allocate(mask_vupdateU(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_vupdateU)'

   allocate(mask_vupdateV(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_vupdateV)'

   allocate(mask_ufluxV(_IRANGE_HALO_-1,_JRANGE_HALO_),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_ufluxV)'

   allocate(mask_vfluxU(_IRANGE_HALO_,_JRANGE_HALO_-1),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (mask_vfluxU)'

   allocate(flux(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (flux)'

   allocate(Di(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Di)'

   allocate(adv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (adv)'
#endif

   mask_updateH  = (az.eq.1)
   mask_uupdateU = (au.eq.1)
   mask_uupdateV = (av.eq.1 .or. av.eq.2)
   mask_vupdateU = (au.eq.1 .or. au.eq.2)
   mask_vupdateV = (av.eq.1)

   mask_ufluxV = ( ax(_IRANGE_HALO_-1,_JRANGE_HALO_).eq.1 )
   mask_vfluxU = ( ax(_IRANGE_HALO_,_JRANGE_HALO_-1).eq.1 )

   adv_gridH%mask_uflux    => mask_vupdateU(_IRANGE_HALO_-1,_JRANGE_HALO_)
   adv_gridH%mask_uupdate  => mask_updateH
   adv_gridH%mask_vflux    => mask_uupdateV(_IRANGE_HALO_,_JRANGE_HALO_-1)
   adv_gridH%mask_vupdate  => mask_updateH
   adv_gridH%mask_finalise => mask_updateH
   adv_gridH%az            => az

   adv_gridU%mask_uflux    => mask_updateH(1+_IRANGE_HALO_,_JRANGE_HALO_)
   adv_gridU%mask_uupdate  => mask_uupdateU
   adv_gridU%mask_vflux    => mask_vfluxU
   adv_gridU%mask_vupdate  => mask_vupdateU ! now also includes y-advection of u along W/E open bdys
   adv_gridU%mask_finalise => mask_vupdateU
   adv_gridU%az            => au

   adv_gridV%mask_uflux    => mask_ufluxV
   adv_gridV%mask_uupdate  => mask_uupdateV ! now also includes x-advection of v along N/S open bdys
   adv_gridV%mask_vflux    => mask_updateH(_IRANGE_HALO_,1+_JRANGE_HALO_)
   adv_gridV%mask_vupdate  => mask_vupdateV
   adv_gridV%mask_finalise => mask_uupdateV
   adv_gridV%az            => av

#if defined(SPHERICAL) || defined(CURVILINEAR)
   adv_gridH%dxu   => dxu(_IRANGE_HALO_-1,_JRANGE_HALO_)
   adv_gridH%dyu   => dyu(_IRANGE_HALO_-1,_JRANGE_HALO_)
   adv_gridH%dxv   => dxv(_IRANGE_HALO_,_JRANGE_HALO_-1)
   adv_gridH%dyv   => dyv(_IRANGE_HALO_,_JRANGE_HALO_-1)
   adv_gridH%arcd1 => arcd1

   adv_gridU%dxu   => dxc(1+_IRANGE_HALO_,_JRANGE_HALO_)
   adv_gridU%dyu   => dyc(1+_IRANGE_HALO_,_JRANGE_HALO_)
   adv_gridU%dxv   => dxx(_IRANGE_HALO_,_JRANGE_HALO_-1)
   adv_gridU%dyv   => dyx(_IRANGE_HALO_,_JRANGE_HALO_-1)
   adv_gridU%arcd1 => arud1

   adv_gridV%dxu   => dxx(_IRANGE_HALO_-1,_JRANGE_HALO_)
   adv_gridV%dyu   => dyx(_IRANGE_HALO_-1,_JRANGE_HALO_)
   adv_gridV%dxv   => dxc(_IRANGE_HALO_,1+_JRANGE_HALO_)
   adv_gridV%dyv   => dyc(_IRANGE_HALO_,1+_JRANGE_HALO_)
   adv_gridV%arcd1 => arvd1
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_advection()'
   write(debug,*)
#endif
   return
   end subroutine init_advection
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_advection - 2D advection schemes \label{sec-do-advection}
!
! !INTERFACE:
   subroutine do_advection(dt,f,U,V,DU,DV,Do,Dn,scheme,split,AH,tag, &
                           Dires,advres)
!
! !DESCRIPTION:
!
! Here, advection terms for all two-dimensional state variables are
! calculated by means of a finite-volume approach (an exception
! is the possibility to directly calculate the momentum advection
! by a one-step two-dimensional upstream scheme,
! see {\tt uv\_advect}) and the advection step is carried out
! as a fractional advection time step. Those 2D variables may be defined on
! T-, U- and V-points. Inside this advection
! routine, it does not matter, where the advected variable is located
! on the grid. All finite volume fluxes and geometric coefficients
! need to be calculated before {\tt do\_advection\_2d} is called.
!
! With the compiler option {\tt SLICE\_MODEL}, the advection in
! $y$-direction is not executed.
!
!
! !USES:
   use halo_zones, only: update_2d_halo,wait_halo,D_TAG,H_TAG,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_ADV,TIM_ADVH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                               :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,Do,Dn,DU,DV
   integer,intent(in)                                :: split,scheme,tag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)        :: f
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: Dires,advres
!
! !LOCAL VARIABLES:
   type(t_adv_grid),pointer :: adv_grid
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection() # ',Ncall
#endif
   call tic(TIM_ADV)

   select case (tag)
      case(H_TAG,D_TAG)
         adv_grid => adv_gridH
      case(U_TAG)
         adv_grid => adv_gridU
      case(V_TAG)
         adv_grid => adv_gridV
      case default
         stop 'do_advection: tag is invalid'
   end select

   Di = Do
   adv = _ZERO_

   select case (split)

      case(NOSPLIT)

         select case (scheme)

            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                _ONE_,scheme,AH,                           &
                                adv_grid%mask_uflux,adv_grid%mask_uupdate, &
#ifdef SLICE_MODEL
                                nosplit_finalise=.true.,                   &
                                mask_finalise=adv_grid%mask_finalise)
#else
                                nosplit_finalise=.false.)
#endif
#ifndef SLICE_MODEL
               call adv_v_split(dt,f,Di,adv,V,Do,DV,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                _ONE_,scheme,AH,                           &
                                adv_grid%mask_vflux,adv_grid%mask_vupdate, &
                                nosplit_finalise=.true.,                   &
                                mask_finalise=adv_grid%mask_finalise)
#endif

            case(UPSTREAM_2DH)

               call adv_upstream_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                     adv_grid%dxv,adv_grid%dyu,   &
                                     adv_grid%dxu,adv_grid%dyv,   &
                                     adv_grid%arcd1,              &
#endif
                                     adv_grid%az,AH)

            case(FCT)

               call adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxv,adv_grid%dyu,   &
                                adv_grid%dxu,adv_grid%dyv,   &
                                adv_grid%arcd1,              &
#endif
                                adv_grid%az,AH)

            case default

               stop 'do_advection: scheme is invalid'

         end select

      case(FULLSPLIT)

         select case (scheme)

            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                _ONE_,scheme,AH,                           &
                                adv_grid%mask_uflux,adv_grid%mask_uupdate)
#ifndef SLICE_MODEL
#ifdef GETM_PARALLEL
               if (scheme.ne.UPSTREAM .and. tag.eq.V_TAG) then
!                 we need to update f(imin:imax,jmax+HALO)
                  call tic(TIM_ADVH)
                  call update_2d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
                  call toc(TIM_ADVH)
               end if
#endif
               call adv_v_split(dt,f,Di,adv,V,Do,DV,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                _ONE_,scheme,AH,                           &
                                adv_grid%mask_vflux,adv_grid%mask_vupdate)
#endif

            case((UPSTREAM_2DH),(FCT))

               stop 'do_advection: scheme not valid for split'

            case default

               stop 'do_advection: scheme is invalid'

         end select

      case(HALFSPLIT)

         select case (scheme)

            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                _HALF_,scheme,AH,                          &
                                adv_grid%mask_uflux,adv_grid%mask_uupdate)
#ifndef SLICE_MODEL
#ifdef GETM_PARALLEL
               if (scheme.ne.UPSTREAM .and. tag.eq.V_TAG) then
!                 we need to update f(imin:imax,jmax+HALO)
                  call tic(TIM_ADVH)
                  call update_2d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
                  call toc(TIM_ADVH)
               end if
#endif
               call adv_v_split(dt,f,Di,adv,V,Do,DV,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxv,adv_grid%dyv,adv_grid%arcd1,  &
#endif
                                _ONE_,scheme,AH,                           &
                                adv_grid%mask_vflux,adv_grid%mask_vupdate)
#endif
#ifdef GETM_PARALLEL
               if (scheme .eq. UPSTREAM) then
                  if (tag .eq. U_TAG) then
!                    we need to update f(imax+1,jmin:jmax)
!                    KK-TODO: if external DU was halo-updated this halo-update is not necessary
                     call tic(TIM_ADVH)
                     call update_2d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,H_TAG)
                     call wait_halo(H_TAG)
                     call toc(TIM_ADVH)
                  end if
               else
!                 we need to update f(imin-HALO:imin-1,jmin:jmax)
!                 we need to update f(imax+1:imax+HALO,jmin:jmax)
                  call tic(TIM_ADVH)
                  call update_2d_halo(f,f,adv_grid%az,imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
                  call toc(TIM_ADVH)
               end if
#endif
               call adv_u_split(dt,f,Di,adv,U,Do,DU,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                adv_grid%dxu,adv_grid%dyu,adv_grid%arcd1,  &
#endif
                                _HALF_,scheme,AH,                          &
                                adv_grid%mask_uflux,adv_grid%mask_uupdate)

            case((UPSTREAM_2DH),(FCT))

               stop 'do_advection: scheme not valid for split'

            case default

               stop 'do_advection: scheme is invalid'

         end select

      case default

         stop 'do_advection: split is invalid'

   end select

   if (present(Dires)) Dires = Di
   if (present(advres)) advres = adv

   call toc(TIM_ADV)
#ifdef DEBUG
   write(debug,*) 'Leaving do_advection()'
   write(debug,*)
#endif
   return
   end subroutine do_advection
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  print_adv_settings
!
! !INTERFACE:
   subroutine print_adv_settings(split,scheme,AH)
!
! !DESCRIPTION:
!
! Checks and prints out settings for 2D advection.
!
! !USES
   IMPLICIT NONE

! !INPUT PARAMETERS:
   integer,intent(in)  :: split,scheme
   REALTYPE,intent(in) :: AH
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'print_adv_settings() # ',Ncall
#endif

   select case (split)
      case((NOSPLIT),(FULLSPLIT),(HALFSPLIT))
      case default
         FATAL 'adv_split=',split,' is invalid'
         stop
   end select

   select case (scheme)
      case((UPSTREAM),(UPSTREAM_2DH),(P2),(SUPERBEE),(MUSCL),(P2_PDM),(FCT))
      case default
         FATAL 'adv_scheme=',scheme,' is invalid'
         stop
   end select

   select case (split)
      case((FULLSPLIT),(HALFSPLIT))
         select case (scheme)
            case((UPSTREAM_2DH),(FCT))
               FATAL 'adv_scheme=',scheme,' not valid for adv_split=',split
               stop
         end select
   end select

   LEVEL3 trim(adv_splits(split))
   LEVEL3 ' ',trim(adv_schemes(scheme))

   if (AH .gt. _ZERO_) then
      LEVEL3 ' with AH=',AH
   else
      LEVEL3 ' without diffusion'
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving print_adv_settings()'
   write(debug,*)
#endif
   return
   end subroutine print_adv_settings
!EOC
!-----------------------------------------------------------------------

   end module advection

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
