#include "cppdefs.h"
!!-----------------------------------------------------------------------
!!BOI
!!
!! !TITLE: Documentation of getm
!!
!! !AUTHORS: Hans Burchard and Karsten Bolding
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
#ifdef _GETM_ESMF_EXEC_
   use getm_esmf, only: do_getm_esmf
#else
#ifdef _GETM_OASIS_
   use getm_oasis, only: do_getm_oasis
#endif
#endif
   use initialise, only: init_model,runtype,dryrun
   use time, only: simtime
   use domain, only: calc_points
   use m2d, only: mem2d
   use getm_timers, only: write_getm_timers
#ifndef NO_3D
   use m3d, only: mem3d
#endif
   use integration
#ifdef GETM_PARALLEL
   use halo_mpi, only: all_2d_exchange, all_3d_exchange
#endif
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! ! LOCAL VARIABLES
   character(len=8)          :: datestr
   character(len=10)         :: timestr
   real                      :: t1=-1,t2=-1,secs
   integer                   :: ierr
!EOP
!-----------------------------------------------------------------------
!BOC
   call cmdline

#ifdef FORTRAN95
   call CPU_Time(t1)
#endif
   call Date_And_Time(datestr,timestr)

#ifdef _GETM_ESMF_EXEC_
   call do_getm_esmf()
#else
#ifdef _GETM_OASIS_
   call do_getm_oasis()
#else
   call init_model(datestr,timestr)
   if ( .not. dryrun ) then
      call time_loop(runtype)
   end if
   call clean_up(dryrun,runtype,MaxN)
#endif
#endif

#ifdef FORTRAN95
   call CPU_Time(t2)
#endif
   call Date_And_Time(datestr,timestr)
   secs = t2-t1
   STDERR LINE
   if( dryrun ) then
      LEVEL1 'getm just did a dry run'
      LEVEL1 'Number of calc-points: ',Calc_Points
      LEVEL1 'Space requirements (global 2D and 3D arrays):'
      LEVEL2 '2D: ',mem2d/1024,' kbytes'
#ifndef NO_3D
      if(runtype .ge. 2) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
#endif
   else
      LEVEL1 'getm: Completed on ',datestr,' ',timestr
      LEVEL1 'Memory used (global 2D and 3D arrays):'
      LEVEL2 '2D: ',mem2d/1024,' kbytes'
#ifndef NO_3D
      if(runtype .ge. 2) then
         LEVEL2 '3D: ',mem3d/1024,' kbytes'
      end if
#endif
      LEVEL1 'Total CPU-time was:    ',secs,' seconds'
      LEVEL1 'Number of time steps:  ',MaxN-MinN+1
      LEVEL1 'Number of calc-points: ',Calc_Points
      if(MaxN-MinN+1 .gt. _ZERO_) then
         LEVEL1 'CPU-time/calc-point:   ',secs/(MaxN-MinN+1)/Calc_Points,' seconds'
         LEVEL1 'Sim-time/CPU-time:     ',simtime/secs
      end if
#ifndef NO_TIMERS
      STDERR LINE
      call write_getm_timers
#endif
   endif
   STDERR LINE
#ifdef GETM_PARALLEL
   LEVEL1 "Communication with other sub-domains:"
   LEVEL2 "2D data exchange: ",all_2d_exchange/(1024*1024)," MB"
   LEVEL2 "3D data exchange: ",all_3d_exchange/(1024*1024)," MB"
   STDERR LINE
#endif
   LEVEL1 'Copyright (C) Karsten Bolding and Hans Burchard.'
   LEVEL1 'under the General Public License (GPL) - http://www.gnu.org '
   STDERR LINE

   !call compilation_options

   end program getm

!EOC

!-----------------------------------------------------------------------
   subroutine cmdline
   use initialise, only: dryrun
   IMPLICIT NONE
   character(len=64)    :: arg
   integer              :: i

   if (command_argument_count() .eq. 0) return

   do i = 1, command_argument_count()
      call get_command_argument(i, arg)

      select case (arg)
      case ('-v', '--version')
         LEVEL0 'GETM www.getm.eu'
         call print_version()
         stop
      case ('-c', '--compile')
         LEVEL0 'GETM www.getm.eu'
         call print_version()
         call compilation_options()
         LEVEL0
         stop
      case ('-h', '--help')
         call print_help()
         stop
!KB      case ('--dryrun')
!KB         dryrun=.true.
      case default
         LEVEL0
         LEVEL0 'Unrecognized command-line option: ', arg
         LEVEL0
         call print_help()
         stop
      end select
   end do
   return
   end

!-----------------------------------------------------------------------
   subroutine print_help()
     character(len=255) :: cmd
     call get_command_argument(0, cmd)

     print '(a)', ''
     print '(a,a,a)', 'usage: ',trim(cmd),' [OPTIONS]'
     print '(a)', ''
     print '(a)', 'Without any options, getm will continue execution.'
     print '(a)', ''
     print '(a)', 'cmdline options:'
     print '(a)', ''
     print '(a)', '  -v, --version     print version information and exit'
     print '(a)', '  -c, --compile     print compilation options'
     print '(a)', '  -h, --help        print usage information and exit'
     print '(a)', ''
     print '(a)', 'visit getm.eu for further info'
     print '(a)', 'consider subscribing to getm-users@googlegroups.com'
     print '(a)', ''
  end subroutine print_help

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------

