   program test_parser
!
! !DESCRIPTION:
!  program to test the parser routine in parser.F90
!  To execute:
!  make test_parser
!
! !USES:
   use parser
   implicit none
! !LOCAL VARIABLES
   integer :: imin=1,imax=100
   integer :: jmin=100,jmax=200
   integer :: kmin=1,kmax=15
   integer :: nmin=1,nmax=10000
!EOP
!-----------------------------------------------------------------------
!BOC
   call parse_output_list('out.dat',imin,imax,jmin,jmax,kmin,kmax,nmin,nmax, &
                         output_list)
   call print_output_list(output_list)

   end program test_parser
!EOC

