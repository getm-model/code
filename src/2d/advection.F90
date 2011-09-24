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
   public init_advection,do_advection
   public adv_u_split,adv_v_split,adv_upstream_2dh,adv_fct_2dh

!  Note (KK): flux is used from the advection routines
!             Di and adv are used only in do_advection
#ifdef STATIC
   REALTYPE,public,dimension(E2DFIELD)        :: flux
   REALTYPE,dimension(E2DFIELD)               :: Di,adv
#else
   REALTYPE,public,dimension(:,:),allocatable :: flux
   REALTYPE,dimension(:,:),allocatable        :: Di,adv
#endif
   integer,public,parameter           :: NOSPLIT=0,FULLSPLIT=1,HALFSPLIT=2
   character(len=64),public,parameter :: adv_splits(0:2) = &
     (/"no split: one 2D uv step",                         &
       "full step splitting: u + v",                       &
       "half step splitting: u/2 + v + u/2"/)
   integer,public,parameter           :: UPSTREAM=1,UPSTREAM_2DH=2,P2=3
   integer,public,parameter           :: SUPERBEE=4,MUSCL=5,P2_PDM=6,FCT=7
   character(len=64),public,parameter :: adv_schemes(7) =  &
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
      subroutine adv_u_split(dt,f,Di,adv,U,Do,DU,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxu,dyu,arcd1,                             &
#endif
                             au,splitfac,scheme,az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,Do,DU
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxu,dyu,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: au,az
         integer,intent(in)                         :: scheme
         logical,intent(in),optional                :: onestep_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_u_split

      subroutine adv_v_split(dt,f,Di,adv,V,Do,DV,                       &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyv,arcd1,                             &
#endif
                             av,splitfac,scheme,az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: V,Do,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: av,az
         integer,intent(in)                         :: scheme
         logical,intent(in),optional                :: onestep_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_v_split

      subroutine adv_upstream_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                  dxv,dyu,dxu,dyv,arcd1,       &
#endif
                                  az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyu,dxu,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         logical,intent(in),optional                :: onestep_finalise
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
      end subroutine adv_upstream_2dh

      subroutine adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                             dxv,dyu,dxu,dyv,arcd1,       &
#endif
                             az,AH,onestep_finalise)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
         REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyu,dxu,dyv,arcd1
#endif
         integer,dimension(E2DFIELD),intent(in)     :: az
         logical,intent(in),optional                :: onestep_finalise
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
! !DESCRIPTION:
!
! Here, memory for some variables is allocated, which are then initialised to
! zero.
!
! !USES
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
   allocate(flux(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (flux)'

   allocate(Di(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Di)'

   allocate(adv(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (adv)'
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
   subroutine do_advection(dt,f,U,V,DU,DV,Do,Dn,                  &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                           dxu,dxv,dyu,dyv,arcd1,                 &
#endif
                           az,au,av,scheme,split,AH,Dires,advres)
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
   use halo_zones, only: update_2d_halo,wait_halo,D_TAG
   use getm_timers, only: tic,toc,TIM_ADV,TIM_ADVH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                               :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(E2DFIELD),intent(in)           :: dxu,dxv,dyu,dyv,arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)            :: az,au,av
   integer,intent(in)                                :: split,scheme
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)        :: f
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

   select case (split)

      case(NOSPLIT)

         select case (scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,      &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,            &
#endif
                                au,_ONE_,scheme,az,AH,    &
#ifdef SLICE_MODEL
                                onestep_finalise=.true.)
#else
                                onestep_finalise=.false.)
#endif

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,Di,adv,V,Do,DV,     &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,           &
#endif
                                av,_ONE_,scheme,az,AH,   &
                                onestep_finalise=.true.)
#endif

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
               stop 'do_advection: scheme=',scheme,' is invalid'
         end select

      case(FULLSPLIT)

         select case (scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,   &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,         &
#endif
                                au,_ONE_,scheme,az,AH)

#ifndef SLICE_MODEL
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
               call adv_v_split(dt,f,Di,adv,V,Do,DV,   &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,         &
#endif
                                av,_ONE_,scheme,az,AH)
#endif

            case((UPSTREAM_2DH),(FCT))
               stop 'do_advection: scheme=',scheme,' not valid for split=',split
            case default
               stop 'do_advection: scheme=',scheme,' is invalid'
         end select

      case(HALFSPLIT)

         select case (scheme)
            case((UPSTREAM),(P2),(SUPERBEE),(MUSCL),(P2_PDM))

               call adv_u_split(dt,f,Di,adv,U,Do,DU,    &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,          &
#endif
                                au,_HALF_,scheme,az,AH)

#ifndef SLICE_MODEL
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
               call adv_v_split(dt,f,Di,adv,V,Do,DV,   &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxv,dyv,arcd1,         &
#endif
                                av,_ONE_,scheme,az,AH)
#endif

               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
               call adv_u_split(dt,f,Di,adv,U,Do,DU,    &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                dxu,dyu,arcd1,          &
#endif
                                au,_HALF_,scheme,az,AH)

            case((UPSTREAM_2DH),(FCT))

               stop 'do_advection: scheme=',scheme,' not valid for split=',split

            case default

               stop 'do_advection: scheme=',scheme,' is invalid'

         end select

      case default

         stop 'do_advection: split=',split,' is invalid'

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
   integer,intent(in)   :: split,scheme
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
