!$Id: test_varinfo.F90,v 1.1 2009-10-06 11:42:49 kb Exp $
   program test_varinfo
!
! !DESCRIPTION:
!  Test program for testing routines in variable_info.F90
!  To execute:
!  make test_varinfo
!
! !USES:
   use variable_info, only: init_var_info, print_var_info
   IMPLICIT NONE
!
! ! LOCAL VARIABLES
!EOP
!-----------------------------------------------------------------------
!BOC
   call init_var_info()
   call print_var_info()

   end program test_varinfo
!EOC
