#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_field - read a 3D field from a file.
!
! !INTERFACE:
   subroutine get_field(fname,var,n,f)
!
! !DESCRIPTION:
!  Reads from file - fname - the variable var into f.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname,var
   integer, intent(in)                 :: n
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer, parameter        :: fmt=NETCDF
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'get_field() # ',ncall
#endif
   select case (fmt)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'Should get an ASCII field'
         stop 'get_field()'
      case (NETCDF)
         call read_field_ncdf(fname,var,n,f)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_field'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_field()'
   write(debug,*)
#endif
   return
   end subroutine get_field
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
