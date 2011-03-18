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
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
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
