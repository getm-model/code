!$Id: read_profile.F90,v 1.1.1.1 2002-05-02 14:01:34 gotm Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: read_profile -
!
! !INTERFACE:
   subroutine read_profile(fn,nmax,zlev,prof,n)
!
! !DESCRIPTION:
!
!  Reads an ASCII file wit a single profile. A very simple format is used
!  first line contains number of elements in the profile - each of the
!  following line lines contains pairs of z-coordinate and variable value.
!  \begin{verbatim}
!   23
!     0.   23.0
!    -2.   22.8
!    :
!    :
!   -50.5   8.6
!  \end{verbatim}
!  The format is free format. Notice the the actual z-coordinate is used
!  with the sea surface being 0. and the z-coordinate positive upwards.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)	:: fn
   integer, intent(in) 	:: nmax
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(inout)	:: zlev(nmax),prof(nmax)
   integer, intent(out) 	:: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: read_profile.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:34  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/08/31 11:22:55  bbh
!  Added n .gt. namx check
!
!  Revision 1.1  2001/08/29 11:11:09  bbh
!  Added read_profile.F90
!
!
! !LOCAL VARIABLES:
   integer 	:: iunit=99,i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   open(unit=iunit,file=trim(fn),err=100,status='unknown')

   read(iunit,*) n

   if (n .gt. nmax) then
       FATAL 'Increase nmax and recompile the calling subroutine'
       stop 'read_profile'
   end if

   do i=n,1,-1
      read(iunit,*) zlev(i),prof(i)
   end do
   close(iunit)

   return
100 FATAL 'Unable to open ',trim(fn)
   stop 'read_profile'
   end subroutine read_profile
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding and Hans Burchard
!-----------------------------------------------------------------------
