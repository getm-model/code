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
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection_3d, do_advection_3d
#ifdef STATIC
   REALTYPE, public                    :: cu(I3DFIELD)
   REALTYPE, public                    :: hi(I3DFIELD)
   REALTYPE, public                    :: hio(I3DFIELD)
#else
   REALTYPE, public, dimension(:,:,:), allocatable       :: hi,hio,cu
#endif
   integer, public, parameter          :: UPSTREAM=1,UPSTREAM_SPLIT=2,P2=3
   integer, public, parameter          :: Superbee=4,MUSCL=5,P2_PDM=6,FCT=7
   REALTYPE, public, parameter         :: one6th=1./6.
   REALTYPE, public, parameter         :: ONE=_ONE_,TWO=2.*_ONE_
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer :: advection_method
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_advection_3d
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

   LEVEL1 'init_advection_3d()'

#ifndef STATIC
   allocate(cu(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (cu)'

   allocate(hi(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hi)'

   allocate(hio(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hio)'
#endif

   cu  = _ZERO_ ; hi  = _ZERO_ ; hio = _ZERO_ 

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
   subroutine do_advection_3d(dt,f,uu,vv,ww,hun,hvn,ho,hn,      &
                             delxu,delxv,delyu,delyv,area_inv,  &
                             az,au,av,hor_adv,ver_adv,adv_split,AH)
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
   use getm_timers, only: tic, toc, TIM_ADVECT3DTOT
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: uu(I3DFIELD)       ! layer-integrated x-transport
   REALTYPE, intent(in) :: vv(I3DFIELD)       ! layer-integrated y-transport
   REALTYPE, intent(in) :: ww(I3DFIELD)       ! grid-related vertical velocity
   REALTYPE, intent(in) :: ho(I3DFIELD)       ! old height of finite volume box
   REALTYPE, intent(in) :: hn(I3DFIELD)       ! new height of finite volume box
   REALTYPE, intent(in) :: hun(I3DFIELD)      ! height of x-interfaces
   REALTYPE, intent(in) :: hvn(I3DFIELD)      ! height of y-interfaces
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
   integer, intent(in)  :: ver_adv            ! selection for vertical scheme
   integer, intent(in)  :: adv_split          ! selection for split mode
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(I3DFIELD)
!
! !LOCAL VARIABLES:
   REALTYPE, parameter       :: a1=0.5*ONE,a2=ONE
   integer         :: k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection_3d() # ',Ncall
#endif
   call tic(TIM_ADVECT3DTOT)

   select case (hor_adv)
      case (UPSTREAM)
         call upstream_adv(dt,f,uu,vv,ww,ho,hn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
      case ((UPSTREAM_SPLIT),(P2),(Superbee),(MUSCL),(P2_PDM),(FCT))
         hi=ho
         select case (adv_split)
            case (0)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a2,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az,& 
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)

#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a2,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az,& 
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
               end if
            case (1)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)

#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
                  call update_3d_halo(f,f,az, &
                                      imin,jmin,imax,jmax,kmax,D_TAG)
                  call wait_halo(D_TAG)

               end if
#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   imin,jmin,imax,jmax,kmax,D_TAG)
               call wait_halo(D_TAG)

            case (2)
               select case (hor_adv)
                  case (UPSTREAM_SPLIT)
                     call upstream_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case (FCT)
                     call fct_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case default
                     FATAL 'For adv_split=2, hor_adv must be 2 (upstream) or 7 (fct)'
               end select
               if (kmax .gt. 1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
               end if
            case default
               FATAL 'Not valid adv_split parameter'
         end select
      case default
         FATAL 'This is not so good - do_advection_3d()'
         stop
   end select

   call toc(TIM_ADVECT3DTOT)
#ifdef DEBUG
   write(debug,*) 'Leaving do_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_advection_3d
!EOC

!-----------------------------------------------------------------------

   end module advection_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
