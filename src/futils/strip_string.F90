#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: strip_string - ....
!
! !INTERFACE:
   subroutine strip_string(a)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  The string a is first stripped for leading and trailing blanks. Then
!  the string is checked for comment characters (presently \! and \#).
!  If a comment character exists all characters to the right are replaced
!  with \<space\>. Then the string is returned. On return the string
!  should be checked with len\_trim() before being used for further
!  processing.
!  The main use of strip_string() is for passing ASCII input files and
!  skipping empty lines and checking for comments before actually
!  using the content of the file.
!
! !USES:
!
! !INOUT PARAMETERS:
   character(len=*), intent(inout)     :: a
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   character(len=255)        :: b
   integer                   :: n,f
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Get rid of leading and trailing blanks
   b = adjustl(a)
   if (len_trim(b) .eq. 0) then
      a = b
      return
   end if
   b = trim(b)

!  check for comment characters - if present overwrite with <space>
   f = scan(b,"#!",back=.false.)
   if (f .gt. 0) then
      do n=f,len(b)
         b(n:n) = ' '
      end do
   end if
   a = b

   return
   end subroutine strip_string
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2015 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
