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
   use halo_zones, only: update_2d_halo,wait_halo,D_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection, do_advection
#ifdef STATIC
   REALTYPE, public                    :: cu(I2DFIELD)
   REALTYPE, public                    :: Di(I2DFIELD)
   REALTYPE, public                    :: Dio(I2DFIELD)
#else
   REALTYPE, public, dimension(:,:), allocatable       :: Di,Dio,cu
#endif
   integer, public, parameter          :: UPSTREAM=1,UPSTREAM_SPLIT=2,P2=3
   integer, public, parameter          :: Superbee=4,MUSCL=5,P2_PDM=6,FCT=7
   character(len=64), public, parameter :: adv_schemes(7) =   &
          (/"2D first-order upstream advection",              &
            "upstream advection (first-order, monotone)",     &
            "P2-PDM advection (third-order, non-monotone)",   &
            "TVD-Superbee advection (second-order, monotone)",&
            "TVD-MUSCL advection (second-order, monotone)",   &
            "TVD-P2-PDM advection (third-order, monotone)",   &
            "2D-FCT advection"/)
   REALTYPE, public, parameter         :: one6th=_ONE_/6
   REALTYPE, public, parameter         :: ONE=_ONE_,TWO=_TWO_
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer :: advection_method
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine adv_u_split(dt,f,U,DU, &
                             delxu,delyu,area_inv,au,splitfac,method,az,AH)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,DU
         REALTYPE,dimension(E2DFIELD),intent(in)    :: delxu,delyu,area_inv
         integer,dimension(E2DFIELD),intent(in)     :: au,az
         integer,intent(in)                         :: method
      end subroutine adv_u_split

      subroutine adv_v_split(dt,f,V,DV, &
                             delxv,delyv,area_inv,av,splitfac,method,az,AH)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,splitfac,AH
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f
         REALTYPE,dimension(E2DFIELD),intent(in)    :: V,DV
         REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyv,area_inv
         integer,dimension(E2DFIELD),intent(in)     :: av,az
         integer,intent(in)                         :: method
      end subroutine adv_v_split

      subroutine adv_upstream(dt,f,U,V,Do,Dn, &
                              delxv,delyu,delxu,delyv,area_inv,az,AH)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn
         REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyu,delxu,delyv
         REALTYPE,dimension(E2DFIELD),intent(in)    :: area_inv
         integer,dimension(E2DFIELD),intent(in)     :: az
      end subroutine adv_upstream

      subroutine adv_upstream_2dh(dt,f,U,V,Do,Dn,DU,DV, &
                                  delxv,delyu,delxu,delyv,area_inv,az,AH)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)                        :: dt,AH
         REALTYPE,dimension(E2DFIELD),intent(inout) :: f
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyu,delxu,delyv
         REALTYPE,dimension(E2DFIELD),intent(in)    :: area_inv
         integer,dimension(E2DFIELD),intent(in)     :: az
      end subroutine adv_upstream_2dh
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
   allocate(cu(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (cu)'

   allocate(Di(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Di)'

   allocate(Dio(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection: Error allocating memory (Dio)'
#endif

   cu  = _ZERO_ ; Di  = _ZERO_ ; Dio = _ZERO_

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
   subroutine do_advection(dt,f,U,V,DU,DV,Do,Dn,             &
                           delxu,delxv,delyu,delyv,area_inv, &
                           az,au,av,hor_adv,adv_split,AH,Dni)
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
   REALTYPE, intent(in) :: U(I2DFIELD)        ! depth-integrated x-transport
   REALTYPE, intent(in) :: V(I2DFIELD)        ! depth-integrated y-transport
   REALTYPE, intent(in) :: Do(I2DFIELD)       ! old depth of water column
   REALTYPE, intent(in) :: Dn(I2DFIELD)       ! new depth of water column
   REALTYPE, intent(in) :: DU(I2DFIELD)      ! depth of x-interfaces
   REALTYPE, intent(in) :: DV(I2DFIELD)      ! depth of y-interfaces
   REALTYPE, intent(in) :: delxu(I2DFIELD)    ! dx centered on u-transport pt.
   REALTYPE, intent(in) :: delxv(I2DFIELD)    ! length of y-interface
   REALTYPE, intent(in) :: delyu(I2DFIELD)    ! length of u-interface
   REALTYPE, intent(in) :: delyv(I2DFIELD)    ! dy centered on v-transport pt.
   REALTYPE, intent(in) :: area_inv(I2DFIELD) ! inverse of horizontal box area
   REALTYPE, intent(in) :: dt                 ! advection time step
   REALTYPE, intent(in) :: AH                 ! constant horizontal diffusivity
   integer, intent(in)  :: az(E2DFIELD)       ! mask for box centre (1: water)
   integer, intent(in)  :: au(E2DFIELD)       ! mask for u-transport (1: water)
   integer, intent(in)  :: av(E2DFIELD)       ! mask for v-transport (1: water)
   integer, intent(in)  :: hor_adv            ! selection for horizontal scheme
   integer, intent(in)  :: adv_split          ! selection for split mode
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(I2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(out), optional   :: Dni(I2DFIELD)
!
! !LOCAL VARIABLES:
   REALTYPE, parameter       :: a1=_HALF_,a2=_ONE_
   integer         :: k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection() # ',Ncall
#endif
   call tic(TIM_ADV)

   select case (hor_adv)
      case (UPSTREAM)
         call adv_upstream(dt,f,U,V,Do,Dn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
         if(present(Dni)) Dni=Dn
      case ((UPSTREAM_SPLIT),(P2),(Superbee),(MUSCL),(P2_PDM),(FCT))
         Di=Do
         select case (adv_split)
            case (0)
               call adv_u_split(dt,f,U,DU,delxu,delyu,area_inv,au,a2,&
                                hor_adv,az,AH)
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,V,DV,delxv,delyv,area_inv,av,a2,&
                                hor_adv,az,AH)
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
#endif

            case (1)
               call adv_u_split(dt,f,U,DU,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

#ifndef SLICE_MODEL
               call adv_v_split(dt,f,V,DV,delxv,delyv,area_inv,av,a2,&
                                hor_adv,az,AH)
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)
#endif

               call adv_u_split(dt,f,U,DU,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call tic(TIM_ADVH)
               call update_2d_halo(f,f,az,imin,jmin,imax,jmax,D_TAG)
               call wait_halo(D_TAG)
               call toc(TIM_ADVH)

            case (2)
               select case (hor_adv)
                  case (UPSTREAM_SPLIT)
                     call adv_upstream_2dh(dt,f,U,V,Do,Dn,DU,DV, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case (FCT)
                     FATAL 'There is a bug in fct_2dh_adv.'
                     !call fct_2dh_adv(dt,f,U,V,Do,Dn,DU,DV, &
                     !          delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case default
                     FATAL 'For adv_split=2, hor_adv must be 2 (upstream) or 7 (fct)'
               end select

            case default
               FATAL 'Not valid adv_split parameter'
         end select
         if(present(Dni)) Dni=Di
      case default
         FATAL 'This is not so good - do_advection()'
         stop
   end select

   call toc(TIM_ADV)
#ifdef DEBUG
   write(debug,*) 'Leaving do_advection()'
   write(debug,*)
#endif
   return
   end subroutine do_advection
!EOC

!-----------------------------------------------------------------------

   end module advection

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
