!$Id: m2d.F90,v 1.2 2003-04-07 12:17:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m2d - depth integrated hydrodynamical model (2D)
!
! !INTERFACE:
   module m2d
!
! !DESCRIPTION:
!  This modules contains declarations for all variables related to 2D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the \emph{domain.F90} module.
!  The module contains public subroutines for initialisation, integration
!  and clean up of the 2D model component.
!  The actual calculation routines are called in integrate\_2d and is linked
!  in from the library lib2d.a.
!
! !USES:
   use time,       only: julianday,secondsofday
   use parameters, only: avmmol
   use domain,     only: imin,imax,jmin,jmax,az,au,av,H,HU,HV,min_depth
   use variables_2d
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG,U_TAG,V_TAG
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   logical	:: have_boundaries
   REALTYPE	:: dtm, z0_const=0.010, Am=-_ONE_
   integer	:: MM=1,residual=-1
   logical 	:: bdy2d=.false.
   integer	:: bdyfmt_2d,bdytype,bdyramp_2d=-1
   character(len=PATH_MAX)	:: bdyfile_2d
   REAL_4B	:: bdy_old(1000)
   REAL_4B	:: bdy_new(1000)
   REAL_4B	:: bdy_data(1000)
   REAL_4B, allocatable :: bdy_times(:)
   integer, parameter	:: comm_method=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: m2d.F90,v $
!  Revision 1.2  2003-04-07 12:17:08  kbk
!  parallel version
!
!  Revision 1.1.1.1  2002/05/02 14:00:41  gotm
!  recovering after CVS crash
!
!  Revision 1.11  2001/10/22 11:55:30  bbh
!  Only call uv_diffusion() if Am > zero
!
!  Revision 1.10  2001/10/22 08:48:30  bbh
!  Am moved from paramters.F90 to m2d.F90
!
!  Revision 1.9  2001/10/11 13:07:03  bbh
!  Support for horizontal diffusion
!
!  Revision 1.8  2001/09/01 17:07:10  bbh
!  Ramping of surface elevation boundaries - via namelist
!
!  Revision 1.7  2001/08/27 11:53:13  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/18 12:55:13  bbh
!  Included masks in calls to update_2d_halo()
!
!  Revision 1.4  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.3  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_2d - initialise 2D relatedstuff.
!
! !INTERFACE:
   subroutine init_2d(runtype,timestep,hotstart)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: runtype
   REALTYPE, intent(in)	:: timestep
   logical, intent(in)	:: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Allocates memiory for 2D related fields.
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer	 	:: rc
   namelist /m2d/ MM,z0_const,Am,residual,bdy2d,bdyfmt_2d,bdyramp_2d,bdyfile_2d
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_2d() # ',Ncall
#endif

   LEVEL1 'init_2d'

!  Read 2D-model specific things from the namelist.
   read(NAMLST,m2d)

   dtm = timestep

!  Allocates memory for the public data members - if not static
   call init_variables_2d(runtype)

#ifdef PARALLEL
   STDERR 'Not calling cfl_check() - PARALLEL'
#else
!KBK   call cfl_check()
#endif

   if (Am .lt. _ZERO_) then
      LEVEL2 'Am is less than zero ---> horizontal diffusion not included'
   end if

   LEVEL2 'Open boundary=',bdy2d
   if (bdy2d) then
      if(hotstart) bdyramp_2d = -1
      LEVEL2 TRIM(bdyfile_2d)
      LEVEL2 'Format=',bdyfmt_2d
   end if

!  Boundary related information
   if (bdy2d) then
!     call have_bdy()
!     call print_bdy('Local Boundary Information')
!kbk     if (have_boundaries) call init_2d_bdy(bdyfmt_2d,bdyfile_2d)
   end if

   call uv_depths()

   where ( -H+min_depth .gt. _ZERO_ )
      z = -H+min_depth
   end where
   zo=z

   where (-HU+min_depth .gt. _ZERO_ )
      zu = -HU+min_depth
   end where
   zub = z0_const ; zub0 = z0_const

   where (-HV+min_depth .gt. _ZERO_ )
      zv = -HV+min_depth
   end where
   zvb = z0_const ; zvb0 = z0_const

   call depth_update()

#ifdef DEBUG
   write(debug,*) 'Leaving init_2d()'
   write(debug,*)
#endif
   return
   end subroutine init_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: integrate_2d - sequence of calls to do 2D model integration
!
! !INTERFACE:
   subroutine integrate_2d(runtype,loop,tausx,tausy,airp)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)		:: runtype,loop
   REALTYPE, intent(in)		:: tausx(E2DFIELD)
   REALTYPE, intent(in)		:: tausy(E2DFIELD)
   REALTYPE, intent(in)		:: airp(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  A wrapper to call all 2D related subroutines in one subroutine.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'integrate_2d() # ',Ncall
#endif

   if (have_boundaries) call update_2d_bdy(loop,bdyramp_2d)

   if (mod(loop-1,MM) .eq. 0) then        ! MacroMicro time step
#ifndef NO_BOTTFRIC
      call bottom_friction(runtype)
#endif
   end if
#ifdef NO_ADVECT
   STDERR 'NO_ADVECT 2D'
#else
#ifndef UV_ADV_DIRECT
   call uv_advect()
   if (Am .gt. _ZERO_) then
      call uv_diffusion(Am) ! Has to be called after uv_advect.
   end if
#endif
#endif
   call momentum(loop,tausx,tausy,airp)
   call sealevel()
   call depth_update()

   if (runtype .gt. 1) then
      Uint=Uint+U
      Vint=Vint+V
   end if
   if(residual .gt. 0 .and. loop .ge. residual) then
      call do_residual(0)
   end if

#ifdef DEBUG
     write(debug,*) 'Leaving integrate_2d()'
     write(debug,*)
#endif
   return
   end subroutine integrate_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_2d - cleanup after 2D run.
!
! !INTERFACE:
   subroutine clean_2d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine cleans up after a 2D integration. Close open files etc.
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
   write(debug,*) 'clean_2d() # ',Ncall
#endif

   call do_residual(1)

#ifdef DEBUG
   write(debug,*) 'Leaving clean_2d()'
   write(debug,*)
#endif
   return
   end subroutine clean_2d
!EOC

!-----------------------------------------------------------------------

   end module m2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
