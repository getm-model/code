!$Id: cleanup.F90,v 1.2 2003-04-07 16:39:16 kbk Exp $
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
   use input, 	only: clean_input
   use meteo, 	only: clean_meteo
   use m2d, 	only: clean_2d
#ifndef NO_3D
   use rivers, 	only: clean_rivers
   use m3d, 	only: clean_3d
#endif
   use output, 	only: clean_output
   use kurt_parallel, 	only: clean_parallel
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
!  Revision 1.2  2003-04-07 16:39:16  kbk
!  parallel support, NO_3D
!
!  Revision 1.1.1.1  2002/05/02 14:01:25  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/09/18 17:48:32  bbh
!  Added algoritm for rivers - getting river data still missing
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
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

   call clean_2d()

#ifndef NO_3D
   call clean_rivers()
   call clean_3d()
#endif

   if( .not. dryrun ) then
      call clean_output()
   end if

   call clean_input()

   call clean_parallel()

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
