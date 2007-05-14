!$Id: adaptive_coordinates.F90,v 1.3.2.1 2007-05-14 11:54:37 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  adaptive vertical coordinates
! \label{sec-adaptive-coordinates}
!
! !INTERFACE:
   subroutine adaptive_coordinates(first)
!
! !DESCRIPTION:
!  For Richard to do
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister and Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

STDERR 'adaptive_coordinates()'

   if (first) then
   end if ! first

#ifdef DEBUG
   write(debug,*) 'Leaving adaptive_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine adaptive_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
