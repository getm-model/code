!$Id: main.F90,v 1.2 2003-03-17 15:02:03 gotm Exp $
#include "cppdefs.h"
!!-----------------------------------------------------------------------
!!BOI
!!
!! !TITLE: Documentation of getm
!!
!! !AUTHORS: Hans Burchard and Karsten Bolding
!!
!! !AFFILIATION: Bolding & Burhard Hydrodynamics
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
   use commhalo, only: myid
   use time,	 only: simtime
   use domain,	 only: calc_points
   use m2d,	 only: mem2d
   use m3d,	 only: mem3d
   use integration
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: main.F90,v $
!  Revision 1.2  2003-03-17 15:02:03  gotm
!  emphasizing (stable version)
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
!
! ! LOCAL VARIABLES
   character(len=8)	:: datestr
   character(len=10)	:: timestr
   real 		:: t1=-1,t2=-1,secs
   integer 		:: ierr
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef FORTRAN95
   call CPU_Time(t1)
#endif
   call Date_And_Time(datestr,timestr)
   STDERR LINE
   STDERR 'getm ver. ',RELEASE,': Started on  ',datestr,' ',timestr
   STDERR LINE

   call compilation_options

   call init_model()
   if ( .not. dryrun ) then
      call time_loop(runtype)
   end if
   call clean_up(dryrun)

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
      if(runtype .gt. 1) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
   else
      LEVEL1 'getm ver. ',RELEASE,': Completed on ',datestr,' ',timestr
      LEVEL1 'Memory used (global 2D and 3D arrays):'
      LEVEL2 '2D: ',mem2d/1024,' kbytes'
      if(runtype .gt. 1) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
      if(secs .gt. _ZERO_) then
         LEVEL1 'Total CPU-time was:    ',secs,' seconds'
      end if
      LEVEL1 'Number of time steps:  ',MaxN-MinN+1
      LEVEL1 'Number of calc-points: ',Calc_Points
      if(secs .gt. _ZERO_) then
         LEVEL1 'CPU-time/calc-point:   ',secs/(MaxN-MinN+1)/Calc_Points,' seconds'
         LEVEL1 'Sim-time/CPU-time:     ',simtime/secs
      end if
   endif
   STDERR LINE
   LEVEL1 'Copyright (C) Bolding & Burchard Hydrodynamics'
   LEVEL1 'under the General Public License (GPL) - http://www.gnu.org '
   STDERR LINE

   if (myid .ge. 0) then
      ierr = 0
   end if

   end program getm

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------


   subroutine compilation_options
   IMPLICIT NONE
!
   STDERR LINE
   STDERR 'Compilation options: (stable version)'
   STDERR LINE
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
#ifdef NO_BOTTFRIC
   LEVEL1 'NO_BOTTFRIC'
#endif
#ifdef NO_ADVECT
   LEVEL1 'NO_ADVECT'
#endif
#ifdef NO_SLR
   LEVEL1 'NO_SLR'
#endif
#ifdef CONST_VISC
   LEVEL1 'CONST_VISC'
#endif
#ifdef PARABOLIC_VISCOSITY
   LEVEL1 'PARABOLIC_VISCOSITY'
#endif
#ifdef MIN_VEL_DEPTH
   LEVEL1 'MIN_VEL_DEPTH'
#endif
#ifdef UV_TVD
   LEVEL1 'UV_TVD'
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

   STDERR LINE

   return
   end
