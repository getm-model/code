!$Id: cleanup.F90,v 1.1 2002-05-02 14:01:25 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  clean_up - a wrapper to perform cleanup
!
! !INTERFACE:
   subroutine clean_up(dryrun)
!
! !DESCRIPTION:
!  Calls the various model components clean up procedures. The called routines
!  should de-allocate memory close open files etc.
!  Some run-time statistics could also be put here.
!
! !USES:
   use output, 	only: clean_output
   use input, 	only: clean_input
   use meteo, 	only: clean_meteo
   use rivers, 	only: clean_rivers
   use m2d, 	only: clean_2d
   use m3d, 	only: clean_3d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)	:: dryrun
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: cleanup.F90,v $
!  Revision 1.1  2002-05-02 14:01:25  gotm
!  Initial revision
!
!  Revision 1.2  2001/09/18 17:48:32  bbh
!  Added algoritm for rivers - getting river data still missing
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_up() # ',Ncall
#endif

   STDERR LINE
   STDERR 'Cleaning up....'
   STDERR LINE

   call clean_meteo()

   call clean_rivers()

   call clean_2d()

   call clean_3d()

   if( .not. dryrun ) then
      call clean_output()
   end if

   call clean_input()

#ifdef DEBUG
   write(debug,*) 'Leaving clean_up()'
   write(debug,*)
#endif

   return
   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
