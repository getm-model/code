!$Id: m3d.F90,v 1.11 2004-01-05 13:23:27 kbk Exp $
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
!  This modules contains declarations for all variables related to 3D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the \emph{domain.F90} module.
!  The module contains public subroutines for initialisation, integration
!  and clean up of the 3D model component.
!  The actual calculation routines are called in integrate\_3d and is linked
!  in from the library lib3d.a.
!
! !USES:
   use parameters, only: avmmol
   use domain, only: maxdepth,vert_cord
   use m2d, only: Am
   use variables_2d, only: D,z,UEx,VEx
#ifndef NO_BAROCLINIC
   use temperature,only: init_temperature, do_temperature
   use salinity,   only: init_salinity, do_salinity
   use eqstate,    only: init_eqstate, do_eqstate
#endif
#ifndef NO_BAROCLINIC
   use suspended_matter, only: init_spm, do_spm
#endif
   use variables_3d
   use advection_3d, only: init_advection_3d
   use bdy_3d, only: init_bdy_3d, do_bdy_3d

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: M=1
   REALTYPE                            :: cord_relax=_ZERO_   
   logical                             :: calc_temp=.true.
   logical                             :: calc_salt=.true.
   logical                             :: calc_spm=.false.
   logical                             :: bdy3d=.false.
   integer                             :: bdyfmt_3d,bdyramp_3d
   character(len=PATH_MAX)             :: bdyfile_3d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: m3d.F90,v $
!  Revision 1.11  2004-01-05 13:23:27  kbk
!  Poor Man's Z-coordinates
!
!  Revision 1.10  2004/01/02 13:54:24  kbk
!  read equation of state info from namelist - Ruiz
!
!  Revision 1.9  2003/12/16 17:02:44  kbk
!  removed TABS - 0. -> _ZERO_
!
!  Revision 1.8  2003/12/16 15:58:54  kbk
!  back ground viscosity and diffusivity (manuel)
!
!  Revision 1.7  2003/09/12 16:23:38  kbk
!  fixed order of use statement - now compiles on with IFC 7.1
!
!  Revision 1.6  2003/08/28 15:13:57  kbk
!  explicit set UEx and VEx to 0
!
!  Revision 1.5  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.4  2003/04/07 16:28:34  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:51  gotm
!  recovering after CVS crash
!
!  Revision 1.30  2001/10/26 09:13:24  bbh
!  Only call slow_diffusion() if Am > 0.
!
!  Revision 1.29  2001/10/23 12:38:58  bbh
!  Forgot cord_relax in one call to coordinates()
!
!  Revision 1.28  2001/10/23 07:53:10  bbh
!  spm -> suspended_matter
!
!  Revision 1.27  2001/10/23 07:05:11  bbh
!  Parabolic viscosity - -DPARABOLIC_VISCOSITY
!
!  Revision 1.26  2001/10/22 11:16:09  bbh
!  Added support for particulate suspended matter - no IO yet
!
!  Revision 1.25  2001/10/22 09:26:41  bbh
!  Added cord_relax
!
!  Revision 1.24  2001/10/22 08:48:30  bbh
!  Am moved from paramters.F90 to m2d.F90
!
!  Revision 1.23  2001/10/12 11:39:20  bbh
!  TVD moved out of ??_momentum_3d.F90 and into uv_advect_3d.F90
!
!  Revision 1.22  2001/10/12 08:49:56  bbh
!  Support for horizontal diffusion
!
!  Revision 1.21  2001/09/19 13:07:00  bbh
!  Moved advection related 3D fields to global allocation
!
!  Revision 1.20  2001/09/03 20:14:03  bbh
!  Fixed reading of vel_{hor,vel,strang
!
!  Revision 1.19  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.18  2001/09/03 13:02:52  bbh
!  Small changes due to NOMADS
!
!  Revision 1.17  2001/08/31 15:39:49  bbh
!  ifdef for CONSTANCE added
!
!  Revision 1.16  2001/08/30 08:56:26  bbh
!  Preparing for 3D boundary conditions
!
!  Revision 1.15  2001/08/29 12:08:05  bbh
!  temp_adv and salt_adv needs to be public
!
!  Revision 1.14  2001/08/29 11:21:46  bbh
!  namelists read in salinity and temperature + initialisation
!
!  Revision 1.13  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.12  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.11  2001/07/26 13:46:06  bbh
!  Included info for salinity and temperature in namelist
!
!  Revision 1.10  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.9  2001/05/25 19:36:36  bbh
!  Added call to eqstate
!
!  Revision 1.8  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.7  2001/05/20 07:51:40  bbh
!  Internal pressure included
!
!  Revision 1.6  2001/05/18 08:17:49  bbh
!  Tidying up initialization of salinity and temperature
!
!  Revision 1.5  2001/05/15 11:47:57  bbh
!  Added update_3d_halo for initial S and T
!
!  Revision 1.4  2001/05/10 11:30:16  bbh
!  Added further support for baroclinicity
!
!  Revision 1.3  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.2  2001/04/24 08:39:21  bbh
!  Included runtype as argument to integrate_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: vel_hor_adv=1,vel_ver_adv=1,vel_strang=0
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_3d - initialise 3D relatedstuff.
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Allocates memiory for 3D related fields.
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: rc
   NAMELIST /m3d/ &
             M,cnpar,cord_relax,                        &
             bdy3d,bdyfmt_3d,bdyramp_3d,bdyfile_3d,     &
             vel_hor_adv,vel_ver_adv,vel_strang,        &
             calc_temp,calc_salt,calc_spm,              &
             avmback,avhback
!
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

! Allocates memory for the public data members - if not static
   call init_variables_3d(runtype)

   LEVEL2 'vel_hor_adv= ',vel_hor_adv
   LEVEL2 'vel_ver_adv= ',vel_ver_adv
   LEVEL2 'vel_strang=  ',vel_strang

   if(vel_hor_adv .gt. 1 .or. vel_ver_adv .gt. 1) then
#ifndef UV_TVD
      STDERR 'To run the model with higher order advection for momentum'
      STDERR 'you need to re-compile the model with the option -DUV_TVD.'
      stop 'init_3d'
#endif
   end if
   dt = M*timestep
   num=avmmol
   nuh=avmmol

!  Needed for interpolation of temperature and salinity
   if (.not. hotstart) then
      call start_macro()
      call coordinates(vert_cord,cord_relax,maxdepth)
   end if

#ifndef NO_BAROCLINIC
   if (runtype .eq. 3 .or. runtype .eq. 4) then
      T = _ZERO_ ; S = _ZERO_ ; rho = _ZERO_
      if(calc_temp) call init_temperature(1)
      if(calc_salt) call init_salinity(1)
      if(calc_spm)  call init_spm(1)
      call init_eqstate()
#ifndef PECS
      call do_eqstate()
#endif
      if (runtype .eq. 3) call internal_pressure()
      if (runtype .eq. 4) then
         call init_advection_3d(2)
      end if
   end if
#endif

   if (bdy3d) call init_bdy_3d()

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
! !IROUTINE: integrate_3d - sequence of calls to do 3D model integration
!
! !INTERFACE:
   subroutine integrate_3d(runtype,n)
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
! A wrapper to call all 3D related subroutines in one subroutine.
!
! !REVISION HISTORY:
!  See log for module
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
   if (bdy3d) call do_bdy_3d(0,T)
#endif

   call coordinates(vert_cord,cord_relax,maxdepth)
#ifndef NO_BOTTFRIC
   if (kmax .gt. 1) then
      call bottom_friction_3d()
   end if
#endif
   SS = _ZERO_
#ifndef NO_BAROCLINIC
   NN = _ZERO_
   if (runtype .eq. 4) call internal_pressure()
#endif
   if (ufirst) then
      call uu_momentum_3d(bdy3d)
      call vv_momentum_3d(bdy3d)
      ufirst=.false.
   else
      call vv_momentum_3d(bdy3d)
      call uu_momentum_3d(bdy3d)
      ufirst=.true.
   end if
   if (kmax .gt. 1) then
      call ww_momentum_3d()
   end if
#ifndef NO_ADVECT
   if (kmax .gt. 1) then
      call uv_advect_3d(vel_hor_adv,vel_ver_adv,vel_strang)
      if (Am .gt. _ZERO_) then
         call uv_diffusion_3d(Am)  ! Must be called after uv_advect_3d
      end if
   end if
#else
   STDERR 'NO_ADVECT 3D'
#endif

#ifdef CONST_VISC
   num = 1.000e-4
   nuh = 1.000e-5
#else
#ifndef NO_BOTTFRIC
   if (kmax .gt. 1) then
      call stresses_3d()
#ifndef PARABOLIC_VISCOSITY
      call ss_nn()
#endif
      call gotm()
   end if
#endif
#endif
#ifndef NO_BAROCLINIC
   if(runtype .eq. 4) then        ! prognostic T and S
      if (calc_temp) call do_temperature(n)
      if (calc_salt) call do_salinity(n)
      if (calc_spm) call do_spm()
#ifndef PECS
      call do_eqstate()
#endif
   end if
#endif

   UEx=_ZERO_ ; VEx=_ZERO_
   if (kmax .gt. 1) then
#ifndef NO_BOTTFRIC
      call slow_bottom_friction()
#endif
#ifndef NO_ADVECT
#ifndef UV_ADV_DIRECT
      call slow_advection()
      if (Am .gt. _ZERO_) then
         call slow_diffusion(Am) ! Has to be called after slow_advection.
      end if
#endif
#endif
   end if

   call slow_terms()
   call stop_macro()

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
! !IROUTINE: clean_3d - cleanup after 3D run.
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
!  This routine cleans up after a 3D integration. Close open files etc.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
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
