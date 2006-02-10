!$Id: stop_macro.F90,v 1.4 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stop_macro() - terminates the macro loop \label{sec-stop-macro}
!
! !INTERFACE:
   subroutine stop_macro
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d, only: Uint,Uinto,Vint,Vinto
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: stop_macro.F90,v $
!  Revision 1.4  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.3  2003-04-23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:56  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
!kbk   integer               :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stop_macro() # ',Ncall
#endif

  Uinto=Uint
  Uint= _ZERO_
  Vinto=Vint
  Vint= _ZERO_

#ifdef DEBUG
   write(debug,*) 'Leaving stop_macro()'
   write(debug,*)
#endif
   return
   end subroutine stop_macro
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
