!$Id: integration.F90,v 1.6 2006-03-17 11:06:32 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  integration - Initialise the time and do the time loop
!
! !INTERFACE:
   module integration
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: MinN=1,MaxN=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: integration.F90,v $
!  Revision 1.6  2006-03-17 11:06:32  kbk
!  cleaner inclusion of SPM module
!
!  Revision 1.5  2005-05-25 10:46:15  kbk
!  introduced progress print out variable - should go in namelist later
!
!  Revision 1.4  2004/03/29 15:35:51  kbk
!  possible to store calculated mean fields
!
!  Revision 1.3  2003/04/23 12:03:46  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 16:39:16  kbk
!  parallel support, NO_3D
!
!  Revision 1.1.1.1  2002/05/02 14:01:25  gotm
!  recovering after CVS crash
!
!  Revision 1.10  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.9  2001/10/17 08:22:56  bbh
!  Resolved conflicts
!
!  Revision 1.8  2001/10/17 08:21:53  bbh
!  Cleaned
!
!  Revision 1.7  2001/09/18 17:48:32  bbh
!  Added algoritm for rivers - getting river data still missing
!
!  Revision 1.6  2001/08/30 08:02:45  bbh
!  Default should be not TEST_METFORCING
!
!  Revision 1.5  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.4  2001/06/04 13:10:25  bbh
!  Cleaning
!
!  Revision 1.3  2001/05/25 19:03:02  bbh
!  New method for - input - all is done via do_input()
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: time_loop - the main loop of getm
!
! !INTERFACE:
   subroutine time_loop(runtype)
!
! !DESCRIPTION:
!  A wrapper that calls meteo\_forcing, integrate\_2d, integrate\_3d and output
!  within a time loop.
!
! !USES:
   use time,     only: update_time,timestep
   use domain,   only: kmax
   use meteo,    only: do_meteo,tausx,tausy,airp
   use m2d,      only: integrate_2d
#ifndef NO_3D
   use m3d,      only: integrate_3d,M
#ifndef NO_BAROCLINIC
   use m3d,      only: T
#endif
   use rivers,   only: do_rivers
#endif
#ifdef SPM
   use suspended_matter, only: spm_calc,do_spm
#endif
   use input,    only: do_input
   use output,   only: do_output,meanout
#ifdef TEST_NESTING
   use nesting,   only: nesting_file
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES
   logical                   :: do_3d
   integer                   :: n
   integer                   :: progress=100
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'time_loop() # ',Ncall
#endif

   STDERR LINE
   LEVEL1 'integrating....'
   STDERR LINE

   do n=MinN,MaxN

      if (progress .gt. 0 .and. mod(n,progress) .eq. 0) LEVEL2 n

#ifndef NO_3D
      do_3d = (runtype .ge. 2 .and. mod(n,M) .eq. 0)
#endif
      call do_input(n)
      if(runtype .le. 2) then
         call do_meteo(n)
#ifndef NO_BAROCLINIC
      else
         call do_meteo(n,T(:,:,kmax))
#endif
      end if
      call integrate_2d(runtype,n,tausx,tausy,airp)
#ifndef NO_3D
      call do_rivers(do_3d)
      if (do_3d) then
         call integrate_3d(runtype,n)

#ifdef SPM
         if (spm_calc) call do_spm()
#endif

#if 0
        call do_biology()
#endif
      end if
#endif

#ifdef TEST_NESTING
      if (mod(n,80) .eq. 0) then
         call nesting_file(WRITING)
      end if
#endif
      call update_time(n)

      if(meanout .ge. 0) then
         call calc_mean_fields(n,meanout)
      end if
      call do_output(runtype,n,timestep)
#ifdef DIAGNOSE
      call diagnose(n,MaxN,runtype)
#endif
   end do

   if (meanout .eq. 0) then
     call calc_mean_fields(n,n)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving time_loop()'
   write(debug,*)
#endif
   return
   end subroutine time_loop
!EOC

!-----------------------------------------------------------------------

   end module integration

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
