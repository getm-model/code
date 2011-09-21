#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 2D advection
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
   use halo_zones, only: update_2d_halo,wait_halo,D_TAG
   IMPLICIT NONE

   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection,do_advection
   public adv_upstream,adv_u_split,adv_v_split,adv_upstream_2dh,adv_fct_2dh

!  Note (KK): cu,Di,adv are used from the advection routines
!             the other fields are provided for uv_advect
#ifdef STATIC
   REALTYPE,public,dimension(E2DFIELD)        :: cu,Di,adv
   REALTYPE,public,dimension(E2DFIELD)        :: Uadv,Vadv,DUadv,DVadv,maskadv,fadv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,public,dimension(E2DFIELD)        :: dxadv,dyadv
#endif
#else
   REALTYPE,public,dimension(:,:),allocatable :: cu,Di,adv
   REALTYPE,public,dimension(:,:),allocatable :: Uadv,Vadv,DUadv,DVadv,maskadv,fadv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,public,dimension(:,:),allocatable :: dxadv,dyadv
#endif
#endif
   integer,public,parameter           :: NOSPLIT=0,FULLSPLIT=1,HALFSPLIT=2
   character(len=64),public,parameter :: adv_splits(0:2) =              &
                                         (/"no split: one 2D step",     &
                                           "splitting: u + v",          &
                                           "splitting: u/2 + v + u/2"/)
   integer,public,parameter           :: UPSTREAM=1,UPSTREAM_2DH=2,P2=3
   integer,public,parameter           :: SUPERBEE=4,MUSCL=5,P2_PDM=6,FCT=7
   character(len=64),public,parameter :: adv_schemes(7) =      &
          (/"upstream advection (first-order, monotone)",      &
            "2DH-upstream advection with forced monotonicity", &
            "P2-PDM advection (third-order, non-monotone)",    &
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
      subroutine adv_u_split(dt,f,Di,adv,U,Do,DU,                           &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxu,dyu,arcd1,                                 &
#endif
                             au,splitfac,adv_scheme,az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,Do,DU
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxu,dyu,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: au,az
         integer,intent(in)                         :: adv_scheme
         logical,intent(in),optional                :: onestep_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_u_split

      subroutine adv_v_split(dt,f,Di,adv,V,Do,DV,                           &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyv,arcd1,                                 &
#endif
                             av,splitfac,adv_scheme,az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: V,Do,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: av,az
         integer,intent(in)                         :: adv_scheme
         logical,intent(in),optional                :: onestep_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_v_split

      subroutine adv_upstream_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                  dxv,dyu,dxu,dyv,arcd1,       &
#endif
                                  az,AH,adv)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyu,dxu,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_upstream_2dh

      subroutine adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyu,dxu,dyv,arcd1,       &
#endif
                             az,AH,adv)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyu,dxu,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_fct_2dh
   end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_advection
!
! !INTERFACE:
   subroutine init_advection()
!
! !DESCRIPTION:
!
! Here, memory for some variables is allocated, which are then initialised to
! zero.
!
! !USES
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection() # ',Ncall
#endif

   LEVEL2 'init_advection()'

#ifndef STATIC
   allocate(cu(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (cu)'

   allocate(Di(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Di)'

   allocate(adv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (adv)'

   allocate(Uadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Uadv)'

   allocate(Vadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Vadv)'

   allocate(DUadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (DUadv)'

   allocate(DVadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (DVadv)'

   allocate(maskadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (maskadv)'

   allocate(fadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (fadv)'

#if defined(SPHERICAL) || defined(CURVILINEAR)
   allocate(dxadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (dxadv)'

   allocate(dyadv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (dyadv)'
#endif
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
   subroutine do_advection(dt,f,U,V,DU,DV,Do,Dn,                          &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                           dxu,dxv,dyu,dyv,arcd1,                         &
#endif
                           az,au,av,adv_scheme,adv_split,AH,Dires,advres)
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
   use getm_timers, only: tic, toc, TIM_ADV, TIM_ADVH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: U(E2DFIELD)        ! depth-integrated x-transport
   REALTYPE, intent(in) :: V(E2DFIELD)        ! depth-integrated y-transport
   REALTYPE, intent(in) :: Do(E2DFIELD)       ! old depth of water column
   REALTYPE, intent(in) :: Dn(E2DFIELD)       ! new depth of water column
   REALTYPE, intent(in) :: DU(E2DFIELD)      ! depth of x-interfaces
   REALTYPE, intent(in) :: DV(E2DFIELD)      ! depth of y-interfaces
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE, intent(in) :: dxu(E2DFIELD)    ! dx centered on u-transport pt.
   REALTYPE, intent(in) :: dxv(E2DFIELD)    ! length of y-interface
   REALTYPE, intent(in) :: dyu(E2DFIELD)    ! length of u-interface
   REALTYPE, intent(in) :: dyv(E2DFIELD)    ! dy centered on v-transport pt.
   REALTYPE, intent(in) :: arcd1(E2DFIELD) ! inverse of horizontal box area
#endif
   REALTYPE, intent(in) :: dt                 ! advection time step
   REALTYPE, intent(in) :: AH                 ! constant horizontal diffusivity
   integer, intent(in)  :: az(E2DFIELD)       ! mask for box centre (1: water)
   integer, intent(in)  :: au(E2DFIELD)       ! mask for u-transport (1: water)
   integer, intent(in)  :: av(E2DFIELD)       ! mask for v-transport (1: water)
   integer, intent(in)  :: adv_scheme         ! selection for horizontal scheme
   integer, intent(in)  :: adv_split          ! selection for split mode
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(E2DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: Dires,advres
!
! !LOCAL VARIABLES:
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

   Di = Do
   adv = _ZERO_

   select case (adv_split)
      case(NOSPLIT)
         select case (adv_scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,         &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,               &
#endif
                                au,_ONE_,adv_scheme,az,AH,   &
#ifdef SLICE_MODEL
                                onestep_finalise=.true.)
#else
                                onestep_finalise=.false.)
#endif

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,Di,adv,V,Do,DV,        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,              &
#endif
                                av,_ONE_,adv_scheme,az,AH,  &
                                onestep_finalise=.true.)
#endif

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

            case(UPSTREAM_2DH)

               call adv_upstream_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                     dxv,dyu,dxu,dyv,arcd1,       &
#endif
                                     az,AH)

            case(FCT)

               call adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyu,dxu,dyv,arcd1,       &
#endif
                                az,AH)

            case default
               FATAL 'Invalid adv_scheme'
               stop
         end select
      case(FULLSPLIT)
         select case (adv_scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,             &
#endif
                                au,_ONE_,adv_scheme,az,AH)

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,Di,adv,V,Do,DV,       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,             &
#endif
                                av,_ONE_,adv_scheme,az,AH)

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
#endif

            case default
               FATAL 'Invalid adv_scheme'
               stop
         end select
      case(HALFSPLIT)
         select case (adv_scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,              &
#endif
                                au,_HALF_,adv_scheme,az,AH)

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,Di,adv,V,Do,DV,       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,             &
#endif
                                av,_ONE_,adv_scheme,az,AH)

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
#endif

               call adv_u_split(dt,f,Di,adv,U,Do,DU,        &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,              &
#endif
                                au,_HALF_,adv_scheme,az,AH)

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

            case default
               FATAL 'Invalid adv_scheme'
               stop
         end select
      case default
         FATAL 'Invalid adv_split'
         stop
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
! !IROUTINE: print_adv_settings
!
! !INTERFACE:
   subroutine print_adv_settings(adv_split,adv_scheme)
!
! !DESCRIPTION:
!
! Checks and prints out settings for advection adv_scheme.
!
! !USES
   IMPLICIT NONE
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

   select case (adv_split)
      case((NOSPLIT),(FULLSPLIT),(HALFSPLIT))
      case default
         FATAL 'Invalid adv_split parameter'
         stop
   end select

   select case (adv_scheme)
      case((UPSTREAM),(UPSTREAM_2DH),(P2),(SUPERBEE),(MUSCL),(P2_PDM),(FCT))
      case default
         FATAL 'Invalid adv_scheme parameter'
         stop
   end select

   select case (adv_split)
      case((FULLSPLIT),(HALFSPLIT))
         select case (adv_scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))
            case((UPSTREAM_2DH),(FCT))
               FATAL 'adv_scheme=',adv_scheme,' not valid for adv_split=',adv_split
               stop
         end select
   end select

   LEVEL3 adv_splits(adv_split)
   LEVEL3 adv_schemes(adv_scheme)

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
