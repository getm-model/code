!$Id: pos.F90,v 1.2 2003-04-23 12:02:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: REALTYPE function pos() - ....
!
! !INTERFACE:
   REALTYPE function pos(a)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: a
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: pos.F90,v $
!  Revision 1.2  2003-04-23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:19  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (a .le. _ZERO_ ) then
      pos = _ZERO_
   else
      pos = a
   end if

   return
   end function pos
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
