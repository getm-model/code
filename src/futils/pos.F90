!$Id: pos.F90,v 1.1.1.1 2002-05-02 14:01:19 gotm Exp $
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
   REALTYPE, intent(in)	:: a
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: pos.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:19  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
!BOC
!kbk   pos=0.5*(a+dabs(a))
!kbk   pos=0.5*(a+abs(a))
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
