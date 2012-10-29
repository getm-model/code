#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m3d - 3D model component
!
! !INTERFACE:
   module m3d
!
! !DESCRIPTION:
!  This module contains declarations for all variables related to 3D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the {\tt domain} module.
!  The module contains public subroutines for initialisation, integration
!  and clean up of the 3D model component.
!  The {\tt m3d} module is initialised in the routine {\tt init\_3d}, see
!  section \ref{sec-init-3d} described on page
!  \pageref{sec-init-3d}.
!  The actual calculation routines are called in {\tt integrate\_3d}
!  (see section \ref{sec-integrate-3d} on page \pageref{sec-integrate-3d}).
!  and are linked in from the library {\tt lib3d.a}.
!  After the simulation, the module is closed in {\tt clean\_3d}, see
!  section \ref{sec-clean-3d} on page \pageref{sec-clean-3d}.
!
! !USES:
   use exceptions
   use parameters, only: avmmol
   use domain, only: openbdy,maxdepth,vert_cord,az
   use m2d, only: uv_advect,uv_diffusion
   use variables_2d, only: z,Uint,Vint,UEx,VEx
#ifndef NO_BAROCLINIC
   use temperature,only: init_temperature, do_temperature, &
            init_temperature_field
   use salinity,   only: init_salinity, do_salinity, init_salinity_field
   use eqstate,    only: init_eqstate, do_eqstate
   use internal_pressure, only: init_internal_pressure, do_internal_pressure
   use internal_pressure, only: ip_method
#endif
   use variables_3d
   use advection, only: NOADV
   use advection_3d, only: init_advection_3d,print_adv_settings_3d,adv_ver_iterations
   use bdy_3d, only: init_bdy_3d, do_bdy_3d
   use bdy_3d, only: bdy3d_tmrlx, bdy3d_tmrlx_ucut, bdy3d_tmrlx_max, bdy3d_tmrlx_min
!  Necessary to use halo_zones because update_3d_halos() have been moved out
!  temperature.F90 and salinity.F90 - should be changed at a later stage
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: M=1
   REALTYPE                            :: cord_relax=_ZERO_
   integer                             :: vel3d_adv_split=0
   integer                             :: vel3d_adv_hor=1
   integer                             :: vel3d_adv_ver=1
   logical                             :: calc_temp=.true.
   logical                             :: calc_salt=.true.
   logical                             :: bdy3d=.false.
   integer                             :: bdyfmt_3d,bdy3d_ramp
   character(len=PATH_MAX)             :: bdyfile_3d
   REALTYPE                            :: ip_fac=_ONE_
   integer                             :: vel_check=0
   REALTYPE                            :: min_vel=-4*_ONE_,max_vel=4*_ONE_
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   logical         :: advect_turbulence=.false.
#ifdef NO_BAROCLINIC
   integer         :: ip_method
#endif
   integer         :: ip_ramp=-1
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_3d - initialise 3D related stuff \label{sec-init-3d}
!
! !INTERFACE:
   subroutine init_3d(runtype,timestep,hotstart)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   REALTYPE, intent(in)                :: timestep
   logical, intent(in)                 :: hotstart
!
!
! !DESCRIPTION:
!  Here, the {\tt m3d} namelist is read from {\tt getm.inp}, and the
!  initialisation of variables is called (see routine {\tt init\_variables}
!  described on page \pageref{sec-init-variables}).
!  Furthermore, a number of consistency checks are made for the choices
!  of the momentum advection schemes. When higher-order advection schemes
!  are chosen for the momentum advection, the compiler option {\tt UV\_TVD}
!  has to be set. Here, the macro time step $\Delta t$ is calculated
!  from the micro time step $\Delta t_m$ and the split factor {\tt M}.
!  Then, in order to have the vertical coordinate system present already here,
!  {\tt coordinates} (see page \pageref{sec-coordinates}) needs to be called,
!  in order to enable proper interpolation of initial values for
!  potential temperature $\theta$ and salinity $S$ for cold starts.
!  Those initial values are afterwards read in via the routines
!  {\tt init\_temperature} (page \pageref{sec-init-temperature}) and
!  {\tt init\_salinity} (page \pageref{sec-init-salinity}).
!  Finally, in order to prepare for the first time step, the momentum advection
!  and internal pressure gradient routines are initialised and the
!  internal pressure gradient routine is called.
!
! !LOCAL VARIABLES:
   integer         :: rc
   NAMELIST /m3d/ &
             M,cnpar,cord_relax,adv_ver_iterations,       &
             bdy3d,bdyfmt_3d,bdy3d_ramp,bdyfile_3d,       &
             bdy3d_tmrlx,bdy3d_tmrlx_ucut,                &
             bdy3d_tmrlx_max,bdy3d_tmrlx_min,             &
             vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver, &
             calc_temp,calc_salt,                         &
             avmback,avhback,advect_turbulence,           &
             ip_method,ip_ramp,                           &
             vel_check,min_vel,max_vel
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_3d() # ',Ncall
#endif

   LEVEL1 'init_3d'
!  Read 3D-model specific things from the namelist.
   read(NAMLST,m3d)
!   rewind(NAMLST)

   LEVEL2 "splitting factor M: ",M

   if (avmback .lt. _ZERO_) then
      LEVEL2 "setting avmback to 0."
      avmback = _ZERO_
   end if
   if (avhback .lt. _ZERO_) then
      LEVEL2 "setting avhback to 0."
      avhback = _ZERO_
   end if

! Allocates memory for the public data members - if not static
   call init_variables_3d(runtype)
   call init_advection_3d()

!  Sanity checks for advection specifications
   LEVEL2 'Advection of horizontal 3D velocities'
#ifdef NO_ADVECT
   if (vel3d_adv_hor .ne. NOADV) then
      LEVEL2 "reset vel3d_adv_hor= ",NOADV," because of"
      LEVEL2 "obsolete NO_ADVECT macro. Note that this"
      LEVEL2 "behaviour will be removed in the future."
      vel3d_adv_hor = NOADV
   end if
   if (vel3d_adv_ver .ne. NOADV) then
      LEVEL2 "reset vel3d_adv_ver= ",NOADV," because of"
      LEVEL2 "obsolete NO_ADVECT macro. Note that this"
      LEVEL2 "behaviour will be removed in the future."
      vel3d_adv_ver = NOADV
   end if
#endif
   call print_adv_settings_3d(vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_)

!  Sanity checks for bdy 3d
   if (.not.openbdy .or. runtype.eq.2) bdy3d=.false.
   if (bdy3d .and. bdy3d_tmrlx) then
      LEVEL2 'bdy3d_tmrlx=.true.'
      LEVEL2 'bdy3d_tmrlx_max=   ',bdy3d_tmrlx_max
      LEVEL2 'bdy3d_tmrlx_min=   ',bdy3d_tmrlx_min
      LEVEL2 'bdy3d_tmrlx_ucut=  ',bdy3d_tmrlx_ucut
      if (bdy3d_tmrlx_min<_ZERO_ .or. bdy3d_tmrlx_min>_ONE_)          &
           call getm_error("init_3d()",                               &
           "bdy3d_tmrlx_min is out of valid range [0:1]")
      if (bdy3d_tmrlx_max<bdy3d_tmrlx_min .or. bdy3d_tmrlx_min>_ONE_) &
           call getm_error("init_3d()",                               &
           "bdy3d_tmrlx_max is out of valid range [bdy3d_tmrlx_max:1]")
      if (bdy3d_tmrlx_ucut<_ZERO_)                                    &
           call getm_error("init_3d()",                               &
           "bdy3d_tmrlx_max is out of valid range [0:inf[")
   end if

   LEVEL2 'ip_ramp=',ip_ramp

   LEVEL2 'vel_check=',vel_check
   if (vel_check .ne. 0) then
      LEVEL3 'doing sanity checks on velocities'
      LEVEL3 'min_vel=',min_vel
      LEVEL3 'max_vel=',max_vel
      if (vel_check .gt. 0) then
         LEVEL3 'out-of-bound values result in termination of program'
      end if
      if (vel_check .lt. 0) then
         LEVEL3 'out-of-bound values result in warnings only'
      end if
   end if

   dt = M*timestep

#ifdef CONSTANT_VISCOSITY
   num=avmback
   nuh=avhback
#else
   num=1.d-15
   nuh=1.d-15

#ifdef TURB_ADV
   if (.not. advect_turbulence) then
      LEVEL2 "reenabled advection of TKE and eps due to"
      LEVEL2 "obsolete TURB_ADV macro. Note that this"
      LEVEL2 "behaviour will be removed in the future!"
      advect_turbulence = .true.
   end if
#endif

   LEVEL2 "advect_turbulence = ",advect_turbulence

#endif

!  Needed for interpolation of temperature and salinity
   if (.not. hotstart) then
      ssen = z
      call start_macro()
      call coordinates(hotstart)
      call hcc_check()
   end if

   if (runtype .eq. 2) then
      calc_temp = .false.
      calc_salt = .false.
   else
#ifndef NO_BAROCLINIC
      T = _ZERO_ ; S = _ZERO_ ; rho = _ZERO_
      if(calc_temp) call init_temperature()
      if(calc_salt) call init_salinity()
#endif
   end if

#ifndef NO_BAROCLINIC
    if (runtype .eq. 3 .or. runtype .eq. 4) then
      call init_eqstate()
#ifndef PECS
      call do_eqstate()
      call ss_nn()
#endif

      if (bdy3d) call init_bdy_3d()
      if (runtype .ge. 3) call init_internal_pressure()
      if (runtype .eq. 3) call do_internal_pressure()
   end if
#endif

   if (vert_cord .eq. _ADAPTIVE_COORDS_) call preadapt_coordinates(preadapt)

#ifdef DEBUG
   write(debug,*) 'Leaving init_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: postinit_3d - re-initialise some 3D after hotstart read.
!
! !INTERFACE:
   subroutine postinit_3d(runtype,timestep,hotstart)
! !USES:
   use domain, only: imin,imax,jmin,jmax, az,au,av
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   REALTYPE, intent(in)                :: timestep
   logical, intent(in)                 :: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine provides possibility to reset/initialize 3D variables to
!  ensure that velocities are correctly set on land cells after read
!  of a hotstart file.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'postinit_3d() # ',Ncall
#endif

   LEVEL1 'postinit_3d'

   if (do_numerical_analyses) then

      allocate(numdis3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'postinit_3d: Error allocating memory (numdis3d)'
      numdis3d = _ZERO_
      allocate(numdis2d(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'postinit_3d: Error allocating memory (numdis2d)'
      numdis2d = _ZERO_

      if (calc_temp) then
         allocate(phymix3d_T(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (phymix3d_T)'
         phymix3d_T = _ZERO_
         allocate(phymix2d_T(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (phymix2d_T)'
         phymix2d_T = _ZERO_
         allocate(nummix3d_T(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (nummix3d_T)'
         nummix3d_T = _ZERO_
         allocate(nummix2d_T(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (nummix2d_T)'
         nummix2d_T = _ZERO_
      end if

      if (calc_salt) then
         allocate(phymix3d_S(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (phymix3d_S)'
         phymix3d_S = _ZERO_
         allocate(phymix2d_S(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (phymix2d_S)'
         phymix2d_S = _ZERO_
         allocate(nummix3d_S(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (nummix3d_S)'
         nummix3d_S = _ZERO_
         allocate(nummix2d_S(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'postinit_3d: Error allocating memory (nummix2d_S)'
         nummix2d_S = _ZERO_
      end if

   end if

! Hotstart fix - see postinit_2d
   if (hotstart) then

      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (au(i,j) .eq. 0) then
               uu(i,j,:)  = _ZERO_
            end if
         end do
      end do
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (av(i,j) .eq. 0) then
               vv(i,j,:)  = _ZERO_
            end if
         end do
      end do
!     These may not be necessary, but we clean up anyway just in case.
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if(az(i,j) .eq. 0) then
               tke(i,j,:) = _ZERO_
               num(i,j,:) = 1.e-15
               nuh(i,j,:) = 1.e-15
#ifndef NO_BAROCLINIC
               S(i,j,:)   = _ZERO_
               T(i,j,:)   = _ZERO_
#endif
            end if
         end do
      end do

      call coordinates(hotstart)

   end if

   return
   end subroutine postinit_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: integrate_3d - calls to do 3D model integration
! \label{sec-integrate-3d}
!
! !INTERFACE:
   subroutine integrate_3d(runtype,n)
   use getm_timers, only: tic, toc, TIM_INTEGR3D
#ifndef NO_BAROCLINIC
   use getm_timers, only: TIM_TEMPH, TIM_SALTH
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype,n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This is a wrapper routine to call all 3D related subroutines.
! The call position for the {\tt coordinates} routine depends on
! the compiler option
! {\tt MUDFLAT}: If it is defined, then the
! call to {\tt coordinates} construction is made such that drying and flooding
! is stable. If {\tt MUDFLAT} is not defined, then the adaptive grids with
! Lagrangian component which are currently under development are supported.
! Both, drying and flooding and
! Lagrangian coordinates does not go together yet.
! The call sequence is as follows:
!
! \vspace{0.5cm}
!
! \begin{tabular}{lll}
! {\tt start\_macro}           & initialising a 3d step & see page
! \pageref{sec-start-macro} \\
! {\tt do\_bdy\_3d}            & boundary conditions for $\theta$ and $S$ & see
! page \pageref{sec-do-bdy-3d} \\
! {\tt coordinates}            & layer heights ({\tt MUTFLAT} defined) & see
! page \pageref{sec-coordinates} \\
! {\tt bottom\_friction\_3d}   & bottom friction & see page
! \pageref{sec-bottom-friction-3d} \\
! {\tt do\_internal\_pressure} & internal pressure gradient & see page
! \pageref{sec-do-internal-pressure} \\
! {\tt uu\_momentum\_3d}       & layer-integrated $u$-velocity & see page
! \pageref{sec-uu-momentum-3d} \\
! {\tt vv\_momentum\_3d}       & layer-integrated $v$-velocity & see page
! \pageref{sec-vv-momentum-3d} \\
! {\tt coordinates}            & layer heights ({\tt MUTFLAT} not defined) & see
! page \pageref{sec-coordinates} \\
! {\tt ww\_momentum\_3d}       & grid-related vertical velocity & see page
! \pageref{sec-ww-momentum-3d} \\
! {\tt uv\_advect\_3d}         & momentum advection & see page
! \pageref{sec-uv-advect-3d} \\
! {\tt uv\_diffusion\_3d}      & momentum diffusion & see page
! \pageref{sec-uv-diffusion-3d} \\
! {\tt stresses\_3d}           & stresses (for GOTM) & see page
! \pageref{sec-stresses-3d} \\
! {\tt ss\_nn}                 & shear and stratification (for GOTM) & see page
! \pageref{sec-ss-nn} \\
! {\tt gotm}                   & interface and call to GOTM & see page
! \pageref{sec-gotm} \\
! {\tt do\_temperature}        & potential temperature equation & see page
! \pageref{sec-do-temperature} \\
! {\tt do\_salinity}           & salinity equation & see page
! \pageref{sec-do-salinity} \\
! {\tt do\_eqstate}            & equation of state & see page
! \pageref{sec-do-eqstate} \\
! {\tt do\_spm}                & suspended matter equation & see page
! \pageref{sec-do-spm} \\
! {\tt do\_getm\_bio}          & call to GOTM-BIO (not yet released) & \\
! {\tt slow\_bottom\_friction} & slow bottom friction & see page
! \pageref{sec-slow-bottom-friction} \\
! {\tt slow\_advection}        & slow advection terms & see page
! \pageref{sec-slow-advection} \\
! {\tt slow\_diffusion}        & slow diffusion terms & see page
! \pageref{sec-slow-diffusion} \\
! {\tt slow\_terms}            & sum of slow terms & see page
! \pageref{sec-slow-terms} \\
! {\tt stop\_macro}            & finishing a 3d step & see page
! \pageref{sec-stop-macro}
! \end{tabular}
!
! \vspace{0.5cm}
!
! Several calls are only executed for certain compiler options. At each
! time step the call sequence for the horizontal momentum equations is
! changed in order to allow for higher order accuracy for the Coriolis
! rotation.
!
! !LOCAL VARIABLES:
  logical, save              :: ufirst=.true.
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'integrate_3d() # ',Ncall
#endif
   call start_macro()
#ifndef NO_BAROCLINIC
#endif
#ifdef MUDFLAT
   call coordinates(.false.)
#endif
#ifndef NO_BOTTFRIC
   if (kmax .gt. 1) then
      call bottom_friction_3d()
   end if
#endif
   call tic(TIM_INTEGR3D)
   SS = _ZERO_
   call toc(TIM_INTEGR3D)
#ifndef NO_BAROCLINIC
   call tic(TIM_INTEGR3D)
   NN = _ZERO_
   call toc(TIM_INTEGR3D)
   if (runtype .eq. 4) call do_internal_pressure()
#endif
   call tic(TIM_INTEGR3D)
   huo=hun
   hvo=hvn
   ip_fac=_ONE_
   if (ip_ramp .gt. 0) ip_fac=min( _ONE_ , n*_ONE_/ip_ramp)
   call toc(TIM_INTEGR3D)
#ifdef STRUCTURE_FRICTION
   call structure_friction_3d
#endif
   if (ufirst) then
      call uu_momentum_3d(n,bdy3d)
      call vv_momentum_3d(n,bdy3d)
      ufirst=.false.
   else
      call vv_momentum_3d(n,bdy3d)
      call uu_momentum_3d(n,bdy3d)
      ufirst=.true.
   end if

!  HERE WE HAVE TO UPDATE SS !!!

#ifndef MUDFLAT
   if (kmax .gt. 1) then
      if (vert_cord .eq. _ADAPTIVE_COORDS_) call ss_nn()
   end if
   call coordinates(.false.)
#endif

   if (kmax .gt. 1) then

      call ww_momentum_3d()

      call uv_advect_3d()
      call uv_diffusion_3d()  ! Must be called after uv_advect_3d

#ifndef NO_BOTTFRIC
      call stresses_3d()
#endif
#ifndef CONSTANT_VISCOSITY
#ifndef PARABOLIC_VISCOSITY
      if (vert_cord .ne. _ADAPTIVE_COORDS_) call ss_nn()
#endif
      call gotm()
      if (advect_turbulence) call tke_eps_advect_3d()
#endif

   end if

#ifndef NO_BAROCLINIC
   if(runtype .eq. 4) then        ! prognostic T and S
      if (calc_temp) call do_temperature(n)
      if (calc_salt) call do_salinity(n)

!     The following is a bit clumpsy and should be changed when do_bdy_3d()
!     operates on individual fields and not as is the case now - on both
!     T and S.
      call tic(TIM_INTEGR3D)
      if (bdy3d) call do_bdy_3d(0,T)
      if (calc_temp) then
         call tic(TIM_TEMPH)
         call update_3d_halo(T,T,az,imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
         call toc(TIM_TEMPH)
         call mirror_bdy_3d(T,D_TAG)
      end if
      if (calc_salt) then
         call tic(TIM_SALTH)
         call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
         call toc(TIM_SALTH)
         call mirror_bdy_3d(S,D_TAG)
      end if
      call toc(TIM_INTEGR3D)

#ifndef PECS
      call do_eqstate()
#endif
!     HERE WE HAVE TO UPDATE NN !!!
   end if
#endif

#ifndef NO_BAROTROPIC
   if (kmax .gt. 1) then
#ifndef NO_BOTTFRIC
      call slow_bottom_friction()
#endif

      call tic(TIM_INTEGR3D)
      call uv_advect(Uint,Vint,Dun,Dvn)
      call uv_diffusion(0,Uint,Vint,Dn,Dun,Dvn) ! Has to be called after uv_advect.
      call toc(TIM_INTEGR3D)

   end if

   call slow_terms()
#endif
   call tic(TIM_INTEGR3D)
   call stop_macro()
   call toc(TIM_INTEGR3D)

#ifdef DEBUG
     write(debug,*) 'Leaving integrate_3d()'
     write(debug,*)
#endif
   return
   end subroutine integrate_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_3d - cleanup after 3D run \label{sec-clean-3d}
!
! !INTERFACE:
   subroutine clean_3d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! Here, a call to the routine {\tt clean\_variables\_3d} which howewer
! does not do anything yet.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_3d() # ',Ncall
#endif

   call clean_variables_3d()

#ifdef DEBUG
     write(debug,*) 'Leaving clean_3d()'
     write(debug,*)
#endif
   return
   end subroutine clean_3d
!EOC

!-----------------------------------------------------------------------

   end module m3d

!-----------------------------------------------------------------------
! Copyright (C) 2000 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
