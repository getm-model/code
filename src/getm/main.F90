!$Id: main.F90,v 1.19 2010-03-03 06:24:22 kb Exp $
#include "cppdefs.h"
!!-----------------------------------------------------------------------
!!BOI
!!
!! !TITLE: Documentation of getm
!!
!! !AUTHORS: Hans Burchard and Karsten Bolding
!!
!! !AFFILIATION: Bolding & Burhard ApS.
!!
!! !DATE:
!!
!! !INTRODUCTION:
!!
!!EOI
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getm - main program
!
! !INTERFACE:
   program getm
!
! !DESCRIPTION:
!
! !USES:
   use initialise, only: init_model,runtype,dryrun
   use time, only: simtime
   use domain, only: calc_points
   use m2d, only: mem2d
   use getm_timers, only: write_getm_timers
#ifndef NO_3D
   use m3d, only: mem3d
#endif
   use integration
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: main.F90,v $
!  Revision 1.19  2010-03-03 06:24:22  kb
!  output info concerning TURB_ADV
!
!  Revision 1.18  2009-09-30 11:28:47  bjb
!  OpenMP threading initial implementation
!
!  Revision 1.17  2009-09-25 12:14:56  kb
!  INCLUDE_HALOS --> SAVE_HALOS
!
!  Revision 1.16  2009-08-21 10:39:00  kb
!  -DINCLUDE_HALOS will include halo-zones when writing/reading NetCDF hotstart files
!
!  Revision 1.15  2009-08-21 08:56:34  bjb
!  Fix name clash on PARALLEL with OpenMP key word
!
!  Revision 1.14  2009-08-18 10:24:46  bjb
!  New getm_timers module
!
!  Revision 1.13  2009-07-28 08:15:42  kb
!  Hydrodynamics --> ApS.
!
!  Revision 1.12  2008-09-16 10:03:24  kb
!  added Holtermanns emergency break algorithm
!
!  Revision 1.11  2008-08-12 08:38:49  kb
!  added NONNEGSALT to compilation_options()
!
!  Revision 1.10  2006-06-02 12:42:20  kbk
!  support for common epoch for hotstart runs
!
!  Revision 1.9  2006-02-06 15:18:21  kbk
!  reverted to v1.7
!
!  Revision 1.7  2005-04-19 15:51:11  kbk
!  notify on use of -DOLD_WRONG_FLUXES
!
!  Revision 1.6  2004/06/15 07:57:49  kbk
!  CONST_VISC --> CONSTANT_VISCOSITY - Ruiz
!
!  Revision 1.5  2003/09/30 09:44:26  kbk
!  hotout=0 -> save hot-files at last time step only
!
!  Revision 1.4  2003/04/23 12:03:46  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:39:16  kbk
!  parallel support, NO_3D
!
!  Revision 1.1.1.1  2002/05/02 14:01:25  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/09/19 14:21:13  bbh
!  Cleaning
!
!  Revision 1.3  2001/09/19 08:34:36  bbh
!  Only calls CPU_time() if -DFORTRAN95
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! ! LOCAL VARIABLES
   character(len=8)          :: datestr
   character(len=10)         :: timestr
   real                      :: t1=-1,t2=-1,secs
   integer                   :: ierr
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef FORTRAN95
   call CPU_Time(t1)
#endif
   call Date_And_Time(datestr,timestr)

   call init_model(datestr,timestr)
   if ( .not. dryrun ) then
      call time_loop(runtype)
   end if
   call clean_up(dryrun,runtype,MaxN)

#ifdef FORTRAN95
   call CPU_Time(t2)
#endif
   call Date_And_Time(datestr,timestr)
   secs = t2-t1
   STDERR LINE
   if( dryrun ) then
      LEVEL1 'getm ver. ',RELEASE,' just did a dry run'
      LEVEL1 'Number of calc-points: ',Calc_Points
      LEVEL1 'Space requirements (global 2D and 3D arrays):'
      LEVEL2 '2D: ',mem2d/1024,' kbytes'
#ifndef NO_3D
      if(runtype .ge. 2) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
#endif
   else
      LEVEL1 'getm ver. ',RELEASE,': Completed on ',datestr,' ',timestr
      LEVEL1 'Memory used (global 2D and 3D arrays):'
      LEVEL2 '2D: ',mem2d/1024,' kbytes'
#ifndef NO_3D
      if(runtype .ge. 2) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
#endif
      if(secs .gt. _ZERO_) then
         LEVEL1 'Total CPU-time was:    ',secs,' seconds'
      end if
      LEVEL1 'Number of time steps:  ',MaxN-MinN+1
      LEVEL1 'Number of calc-points: ',Calc_Points
      if(secs .gt. _ZERO_) then
         LEVEL1 'CPU-time/calc-point:   ',secs/(MaxN-MinN+1)/Calc_Points,' seconds'
         LEVEL1 'Sim-time/CPU-time:     ',simtime/secs
      end if
#ifndef NO_TIMERS
      STDERR LINE
      call write_getm_timers
#endif      
   endif
   STDERR LINE
   LEVEL1 'Copyright (C) Bolding & Burchard ApS.'
   LEVEL1 'under the General Public License (GPL) - http://www.gnu.org '
   STDERR LINE

   call compilation_options

   end program getm

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------


   subroutine compilation_options
   IMPLICIT NONE
!
   STDERR LINE
   STDERR 'Compilation options (unstable version):'
   STDERR LINE
!
#ifdef GETM_PARALLEL
   LEVEL1 'Compiled for parallel execution'
#else
   LEVEL1 'Compiled for serial execution'
#endif

#ifdef GETM_OMP
   LEVEL1 '   with OpenMP thread capability'
#else
   LEVEL1 '   without OpenMP thread capability'
#endif

!
#ifdef NO_3D
   LEVEL1 'NO_3D'
#endif
#ifdef NO_BAROCLINIC
   LEVEL1 'NO_BAROCLINIC'
#endif
!
#ifdef FORTRAN90
   LEVEL1 'Fortran 90 compilation'
#endif
!
#ifdef FORTRAN95
   LEVEL1 'Fortran 95 compilation'
#endif
!
#ifdef PRODUCTION
   LEVEL1 'Production compilation'
#endif
!
#ifdef PROFILING
   LEVEL1 'Profiling is enabled'
#endif
!
#ifdef DEBUG
   LEVEL1 'Debugging enabled'
#endif
!
#ifdef STATIC
   LEVEL1 'Using STATIC memory allocation'
#else
   LEVEL1 'Using DYNAMIC memory allocation'
#endif
!
#ifdef SINGLE
   LEVEL1 'Using single precision'
#else
   LEVEL1 'Using double precision'
#endif
!
! Various tests
#ifdef SPHERICAL
   LEVEL1 'SPHERICAL'
#endif
#ifdef CURVILINEAR
   LEVEL1 'CURVILINEAR'
#endif
#ifdef TURB_ADV
   LEVEL1 'TURB_ADV'
#endif
#ifdef NO_BOTTFRIC
   LEVEL1 'NO_BOTTFRIC'
#endif
#ifdef NO_ADVECT
   LEVEL1 'NO_ADVECT'
#endif
#ifdef NO_SLR
   LEVEL1 'NO_SLR'
#endif
#ifdef CONSTANT_VISCOSITY
   LEVEL1 'CONSTANT_VISCOSITY'
#endif
#ifdef PARABOLIC_VISCOSITY
   LEVEL1 'PARABOLIC_VISCOSITY'
#endif
#ifdef MIN_VEL_DEPTH
   LEVEL1 'MIN_VEL_DEPTH'
#endif
#ifdef NEW_SS
   LEVEL1 'NEW_SS'
#endif
#ifdef UV_TVD
   LEVEL1 'UV_TVD'
#endif
#ifdef NONNEGSALT
   LEVEL1 'NONNEGSALT'
#endif
#ifdef USE_BREAKS
   LEVEL1 'USE_BREAKS'
#endif
#ifdef PRESS_GRAD_Z
   LEVEL1 'PRESS_GRAD_Z'
#endif
#ifdef ITERATE_VERT_ADV
   LEVEL1 'ITERATE_VERT_ADV'
#endif
#ifdef SUBSTR_INI_PRESS
   LEVEL1 'SUBSTR_INI_PRESS'
#endif
#ifdef SONG_WRIGHT
   LEVEL1 'SONG_WRIGHT'
#endif
#ifdef NO_TIMERS
   LEVEL1 'NO_TIMERS'
#endif
#ifdef OLD_WRONG_FLUXES
   LEVEL1 'OLD_WRONG_FLUXES'
#endif
#ifdef _WRITE_HALOS_
   LEVEL1 '_WRITE_HALOS_'
#endif
#ifdef _WRITE_HOT_HALOS_
   LEVEL1 '_WRITE_HOT_HALOS_'
#endif
#ifdef _READ_HOT_HALOS_
   LEVEL1 '_READ_HOT_HALOS_'
#endif

   STDERR LINE

   return
   end
