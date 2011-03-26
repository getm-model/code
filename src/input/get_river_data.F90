#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_river_data - read river data from file.
!
! !INTERFACE:
   subroutine get_river_data(n)
!
! !DESCRIPTION:
!  Reads river data from data file
!
! !USES:
   use ncdf_river, only: get_river_data_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'get_river_data() # ',ncall
#endif

   select case (NETCDF)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'should get ASCII river data'
         stop 'ASCII river data input not coded!'
      case (NETCDF)
         call get_river_data_ncdf(n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_river_data'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_river_data()'
   write(debug,*)
#endif
   return
   end subroutine get_river_data
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
