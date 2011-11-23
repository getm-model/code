#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  getm_timers - use tic-toc principle to time parts of code
!
! !INTERFACE:
   MODULE getm_timers
!
! !DESCRIPTION:
!  This module implements a set of variables and subroutines which make
!  it possible to find cumulative wall-time spent on various parts of
!  GETM. A subroutine can call tic and toc to start and stop timers.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_getm_timers
   public                              :: write_getm_timers,write_getm_timers_active
   public                              :: tic, toc
! !PUBLIC DATA MEMBERS:
! The indices order the output.
! Allow extra "spaces" here and there to let new timers be possible
! without having to renumber the entire list
! Indices into common arrays:                  ! Subroutine to time:
   integer, parameter :: TIM_INITIALIZE  =  1   ! initialize
   integer, parameter :: TIM_BOTTFRICT   =  2   ! 2d bottom_friction
   integer, parameter :: TIM_MOMENTUM    =  4   ! 2d momentum
   integer, parameter :: TIM_MOMENTUMH   =  5   ! 2d momentum - halo part only
   integer, parameter :: TIM_UVDEPTHS    =  6   ! 2d uv_depths
   integer, parameter :: TIM_UVEX        =  7   ! wrapper for hor. adv+diff
   integer, parameter :: TIM_UVADV       =  8   ! 2d uv_advect
   integer, parameter :: TIM_UVADVH      =  9   ! 2d uv_advect - halo part only
   integer, parameter :: TIM_UVDIFF      = 10   ! 2d uv_diffusion
   integer, parameter :: TIM_DPTHUPDATE  = 12   ! 2d depth_update
   integer, parameter :: TIM_SEALEVEL    = 14   ! 2d sealevel
   integer, parameter :: TIM_SEALEVELH   = 15   ! 2d sealevel - halo part only
   integer, parameter :: TIM_GOTM        = 20   ! 3d gotm
   integer, parameter :: TIM_GOTMTURB    = 22   ! 3d gotm - external calls only (do_turbulence)
   integer, parameter :: TIM_GOTMH       = 23   ! 3d gotm - halo part onl
   integer, parameter :: TIM_UVADV3D     = 30   ! 3d uv_advect_3d
   integer, parameter :: TIM_UVADV3DH    = 31   ! 3d uv_advect_3d - halo part only
   integer, parameter :: TIM_UVDIFF3D    = 32   ! 3d uv_diffusion_3d
   integer, parameter :: TIM_VVMOMENTUM  = 34   ! 3d vv_momentum_3d
   integer, parameter :: TIM_VVMOMENTUMH = 35   ! 3d vv_momentum_3d - halo part only
   integer, parameter :: TIM_UUMOMENTUM  = 36   ! 3d uu_momentum_3d
   integer, parameter :: TIM_UUMOMENTUMH = 37   ! 3d uu_momentum_3d - halo part only
   integer, parameter :: TIM_WWMOMENTUM  = 38   ! 3d ww_momentum_3d
   integer, parameter :: TIM_WWMOMENTUMH = 39   ! 3d ww_momentum_3d - halo part only
   integer, parameter :: TIM_SSNN        = 40   ! 3d ss_nn
   integer, parameter :: TIM_BOTTFRICT3D = 42   ! 3d bottom_friction_3d
   integer, parameter :: TIM_STRESSES3D  = 44   ! 3d stresses_3d
   integer, parameter :: TIM_STRESSES3DH = 45   ! 3d stresses_3d - halo part only
   integer, parameter :: TIM_SLOWTERMS   = 46   ! 3d slow_terms
   integer, parameter :: TIM_SLOWBFRICT  = 48   ! 3d slow_bottom_friction
   integer, parameter :: TIM_TEMP        = 52   ! 3d do_temperature
   integer, parameter :: TIM_TEMPH       = 53   ! 3d temperature halo (presently in m3d/do_integrate_3d)
   integer, parameter :: TIM_SALT        = 54   ! 3d do_salinity
   integer, parameter :: TIM_SALTH       = 55   ! 3d salinity halo (presently in m3d/do_integrate_3d)
   integer, parameter :: TIM_COORDS      = 56   ! 3d coordinates
   integer, parameter :: TIM_INTPRESS    = 58   ! 3d do_internal_pressure
   integer, parameter :: TIM_STARTMCR    = 60   ! 3d start_macro
   integer, parameter :: TIM_STOPMCR     = 62   ! 3d stop_macro
   integer, parameter :: TIM_STRCTFRICT  = 64   ! 3d structure_friction_3d
   integer, parameter :: TIM_EQSTATE     = 66   ! 3d do_eqstate
   integer, parameter :: TIM_CALCMEANF   = 68   ! 3d calc_mean_fields
   integer, parameter :: TIM_METEO       = 70   ! do_meteo (could use + halo)
   integer, parameter :: TIM_GETM_BIO    = 72   ! do_getm_bio
   ! These catch compuations in integrate_[23]d, which are not in other timers:
   integer, parameter :: TIM_INTEGR2D    = 80   ! 2d integrate_2d - remaining stuff
   integer, parameter :: TIM_INTEGR3D    = 81   ! 3d integrate_3d - remaining stuff
   ! This is for do_input and do_output
   integer, parameter :: TIM_INPUT       = 90   ! input
   integer, parameter :: TIM_OUTPUT      = 92   ! output
   ! These catch stuff that are *also* measured somewhere else:
   integer, parameter :: TIM_ADV         = 100  ! 2d advection
   integer, parameter :: TIM_ADVH        = 101  ! 2d advection halo parts
   integer, parameter :: TIM_ADV3D       = 102  ! 3d advection
   integer, parameter :: TIM_ADV3DH      = 103  ! 3d advection halo parts
   integer, parameter :: TIM_CHECK3DF    = 104  ! check_3d_fields
   integer, parameter :: TIM_MIXANALYSIS = 105  ! (numerical) mixing analysis
   integer, parameter :: TIM_HALO2D      = 110  ! do halo 2d (initialize comm)
   integer, parameter :: TIM_HALO3D      = 111  ! do halo 3d (initialize comm)
   integer, parameter :: TIM_HALOWAIT    = 112  ! wait_halo (2d+3d both)
   ! This is test timers for temporary coding purposes:
   !  Note: All timers with index 170+ (test_timer_first) are
   !  considered test timers, so dont implement your timers here
   integer, parameter :: TIM_TEST00      = 170
   integer, parameter :: TIM_TEST01      = 171
   integer, parameter :: TIM_TEST02      = 172
   integer, parameter :: TIM_TEST03      = 173
   integer, parameter :: TIM_TEST04      = 174
   integer, parameter :: TIM_TEST05      = 175
   integer, parameter :: TIM_TEST06      = 176
   integer, parameter :: TIM_TEST07      = 177
   integer, parameter :: TIM_TEST08      = 178
   integer, parameter :: TIM_TEST09      = 179
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann

! !LOCAL VARIABLES:
   integer, parameter :: max_timers       = 200
   integer, parameter :: test_timer_first = 170
   LONGINT            :: timercounts(max_timers)
   LONGINT            :: timertics(max_timers), sysclockcalls(max_timers)
   character(len=24)  :: timernames(max_timers)
   LONGINT, save      :: num_clock_calls=0
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_timers - initialise the subroutine timers
!
! !INTERFACE:
   subroutine init_getm_timers()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initializes timers for some subroutines.
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
! !LOCAL VARIABLES:
   integer                   :: i
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_getm_timers() # ',Ncall
#endif

!
! Initialize all timers
   do i=1,max_timers
      timercounts(i)   = 0
      timertics(i)     = 0
      sysclockcalls(i) = 0
      timernames(i)    = ''
   end do

!
! For logging/debugging store name of each subroutine.
! Only the timers with defined names will be shown at end.
!
! It is possible to time less than a single subroutine.
! If you put "timers within timers", then let the first character of
! the "inner" timers start with a blank. This formats well and is
! used when computing total wall time
!
!   The integers are to easily see length of a string compare to max=24.
!                                0         1         2
!                                0123456789012345678901234
   timernames(TIM_INITIALIZE)  = 'initialize'

   timernames(TIM_BOTTFRICT)   = 'bottom_friction'
   timernames(TIM_MOMENTUM)    = 'momentum'
   timernames(TIM_UVDEPTHS)    = 'uv_depths'
   timernames(TIM_UVEX)        = ' sum calc_uvex'
   timernames(TIM_UVADV)       = ' sum uv_advect'
   timernames(TIM_UVDIFF)      = ' sum uv_diffusion'
   timernames(TIM_DPTHUPDATE)  = 'depth_update'
   timernames(TIM_SEALEVEL)    = 'sealevel'

   timernames(TIM_INTEGR2D)    = 'integrate_2d other'
   timernames(TIM_ADV)         = ' sum do_advection'

   timernames(TIM_METEO)       = 'do_meteo'
   timernames(TIM_INPUT)       = 'do_input'
   timernames(TIM_OUTPUT)      = 'do_output'

#ifdef GETM_PARALLEL
   timernames(TIM_MOMENTUMH)   = ' momentum-halo'
   timernames(TIM_UVADVH)      = ' uv_advect-halo'
   timernames(TIM_SEALEVELH)   = ' sealevel-halo'
   timernames(TIM_ADVH)        = ' do_advection-halo'
   timernames(TIM_HALO2D)      = ' sum do_halo_2d'
   timernames(TIM_HALOWAIT)    = ' sum wait_halo'
#endif

#ifndef NO_3D
   timernames(TIM_GOTM)        = 'gotm'
   ! do_turbulence is deeply nested and requires *many* calls to time
   !timernames(TIM_GOTMTURB)    = ' gotm-turbulence'
   timernames(TIM_UVADV3D)     = 'uv_advect_3d'
   timernames(TIM_UVDIFF3D)    = 'uv_diffusion_3d'
   timernames(TIM_VVMOMENTUM)  = 'vv_momentum_3d'
   timernames(TIM_UUMOMENTUM)  = 'uu_momentum_3d'
   timernames(TIM_WWMOMENTUM)  = 'ww_momentum_3d'
   timernames(TIM_SSNN)        = 'ss_nn'
   timernames(TIM_BOTTFRICT3D) = 'bottom_friction_3d'
   timernames(TIM_STRESSES3D)  = 'stresses_3d'
   timernames(TIM_SLOWTERMS)   = 'slow_terms'
   timernames(TIM_SLOWBFRICT)  = 'slow_bottom_friction'
   timernames(TIM_TEMP)        = 'do_temperature'
   timernames(TIM_SALT)        = 'do_salinity'
   timernames(TIM_COORDS)      = 'coordinates'
   timernames(TIM_INTPRESS)    = 'do_internal_pressure'
   timernames(TIM_STARTMCR)    = 'start_macro'
   timernames(TIM_EQSTATE)     = 'eq_state'
   timernames(TIM_CALCMEANF)   = 'calc_mean_fields'

   timernames(TIM_CHECK3DF)    = ' sum check_3d_fields'
   timernames(TIM_ADV3D)       = ' sum do_advection_3d'
   timernames(TIM_INTEGR3D)    = 'integrate_3d other'
   timernames(TIM_MIXANALYSIS) = 'numerical mixing analysis'

! We only really want to display halo-stuff if we compile for parallel:
#ifdef GETM_PARALLEL
   timernames(TIM_GOTMH)       = ' gotm-halo'
   timernames(TIM_UVADV3DH)    = ' uv_advect_3d-halo'
   timernames(TIM_UUMOMENTUMH) = ' uu_momentum_3d-halo'
   timernames(TIM_VVMOMENTUMH) = ' vv_momentum_3d-halo'
   timernames(TIM_WWMOMENTUMH) = ' ww_momentum_3d-halo'
   timernames(TIM_STRESSES3DH) = ' stresses_3d-halo'
   timernames(TIM_TEMPH)       = ' temperature-halo'
   timernames(TIM_SALTH)       = ' salinity-halo'
   timernames(TIM_ADV3DH)      = ' do_advection_3d-halo'
   timernames(TIM_HALO3D)      = ' sum do_halo_3d'
#endif

#ifdef STRUCTURE_FRICTION
   timernames(TIM_STRCTFRICT)  = 'structure_friction_3d'
#endif
#endif

   timernames(TIM_TEST00)  = ' test-00'
   timernames(TIM_TEST01)  = ' test-01'
   timernames(TIM_TEST02)  = ' test-02'
   timernames(TIM_TEST03)  = ' test-03'
   timernames(TIM_TEST04)  = ' test-04'
   timernames(TIM_TEST05)  = ' test-05'
   timernames(TIM_TEST06)  = ' test-06'
   timernames(TIM_TEST07)  = ' test-07'
   timernames(TIM_TEST08)  = ' test-08'
   timernames(TIM_TEST09)  = ' test-09'

#ifdef DEBUG
   write(debug,*) 'Leaving init_timers()'
   write(debug,*)
#endif
   return
   end subroutine init_getm_timers
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: tic - Start wall clock timer for particular subroutine
!
! !INTERFACE:
   subroutine tic(timerindex)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Start/store wall clock for particular timer
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: timerindex
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NO_TIMERS
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tic() # ',Ncall
#endif
!
   CALL SYSTEM_CLOCK(timertics(timerindex))
   num_clock_calls = num_clock_calls+1
   sysclockcalls(timerindex) = sysclockcalls(timerindex)+1
#ifdef DEBUG
   write(debug,*) 'Leaving tic()'
   write(debug,*)
#endif
#endif
   return
   end subroutine tic
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: toc - Stop and record wall clock timer
!
! !INTERFACE:
   subroutine toc(timerindex)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Get wall clock and update used time for particular timer
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: timerindex
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
! !LOCAL VARIABLES:
   LONGINT              :: timertoc
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NO_TIMERS
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'toc() # ',Ncall
#endif

   CALL SYSTEM_CLOCK(timertoc)
   num_clock_calls = num_clock_calls+1
   sysclockcalls(timerindex) = sysclockcalls(timerindex)+1
   if (timertoc .ge. timertics(timerindex)) then
      timercounts(timerindex) =  timercounts(timerindex) + (timertoc-timertics(timerindex))
   else
      LEVEL2 'Warning. System_clock clicks timing ',timernames(timerindex),' rolled past maximum'
      LEVEL2 '  Timings may be slightly incorrect'
   end if

   timertics(timerindex) = 0

#ifdef DEBUG
   write(debug,*) 'Leaving toc()'
   write(debug,*)
#endif
#endif
   return
   end subroutine toc
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_getm_timers - write content of subroutine timers
!
! !INTERFACE:
   subroutine write_getm_timers()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Write contents of all timers
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
! !LOCAL VARIABLES:
   integer                   :: i
   LONGINT                   :: count_rate, count_max, count_dummy
   REALTYPE                  :: thistime,tottime
!EOP
!-------------------------------------------------------------------------
!BOC
#ifndef NO_TIMERS
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'write_getm_timers() # ',Ncall
#endif

!
! Get local tics per second for this machine:
   CALL SYSTEM_CLOCK(count_dummy,count_rate,count_max)
!
! Compute total clocked wall time for later use in output.
   tottime=_ZERO_
   do i=1,max_timers
      if (LGE(timernames(i),'0')) then
         ! "Minor" timers should start with space (ASCII char 0).
         ! Actually, dash or underscore is ok too.
         ! This is not an "inner" timer add to total timed wall time:
         tottime=tottime+_ONE_*timercounts(i)/count_rate
      end if
   end do
!
! Write all timers
   LEVEL1 'GETM timers in seconds (',count_rate,' tics per second)'
   LEVEL1 '                             wall      #sysclocks     % of total'
   LEVEL1 ' Timername                time (s)        calls       clocked time'
   do i=1,max_timers
      if (len_trim(timernames(i)).gt.0) then
         thistime = _ONE_*timercounts(i)/count_rate
!  Note: Test-timers (index test_timer_first+) only written if
!   they are actually used
         if (i.lt.test_timer_first .or. sysclockcalls(i).gt.0) then
            write(stderr,100) timernames(i),thistime,                   &
                 sysclockcalls(i),100*thistime/tottime
100         FORMAT(A26,F10.2,I15,F14.2)
         end if
     end if
   end do
   write(stderr,100) '  TIMERS SUM  ',tottime,num_clock_calls,100*_ONE_

!
! For logging/debugging store name of each subroutine.
! In reality, only these timers should be used.
! Make sure that the first character is *not* blank.

#ifdef DEBUG
   write(debug,*) 'Leaving write_getm_timers()'
   write(debug,*)
#endif
#endif
   return
   end subroutine write_getm_timers
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_getm_timers_active - write which timeres are ON
!
! !INTERFACE:
   subroutine write_getm_timers_active()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Write IDs of active timers.
!  This is mostly useful for debugging and location of where certain
!  timers are set or should be set.
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
! !LOCAL VARIABLES:
   integer                   :: i,nactive
   integer                   :: active_counters(max_timers)
!EOP
!-------------------------------------------------------------------------
!BOC
#ifndef NO_TIMERS
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'write_getm_timers_active() # ',Ncall
#endif

   nactive=0
! Go over each timer and mark it if it is active
   do i=1,max_timers
      if (timertics(i).ne.0) then
         ! Inactive timers have been set to "zero" tic value.
         nactive=nactive+1
         active_counters(nactive)=i
      end if
   end do
   if (nactive>0) then
      LEVEL2 'Active timers:',active_counters(1:nactive)
   else
      LEVEL2 'Active timers: none'
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving write_getm_timers_active()'
   write(debug,*)
#endif
#endif
   return
 end subroutine write_getm_timers_active
!EOC
!-----------------------------------------------------------------------

   end module getm_timers

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Karsten Bolding and Bjarne Buchmann             !
!-----------------------------------------------------------------------
