#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: inquire_file - inquire a file for variable info
!
! !INTERFACE:
   subroutine inquire_file(fn,ncid,varids,varnames)
!
! !DESCRIPTION:
!  Reads varname from a named file - fn - into to field.
!
! !USES:
   use ncdf_get_field, only: inquire_file_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !OUTPUT PARAMETERS:
   integer, intent(inout)              :: ncid
   integer, allocatable, intent(inout) :: varids(:)
   character(len=50), allocatable, intent(inout) :: varnames(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
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
         stop 'get_2d_field()'
      case (NETCDF)
         call inquire_file_ncdf(fn,ncid,varids,varnames)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_2d_field'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving inquire_file()'
   write(debug,*)
#endif
   return
   end subroutine inquire_file
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
