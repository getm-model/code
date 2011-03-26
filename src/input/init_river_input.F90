#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_river_input - initialise 2D boundary data file(s)
!
! !INTERFACE:
   subroutine init_river_input(fn,n)
!
! !DESCRIPTION:
!  Prepares for reading boundary data for the 2D module.
!
! !USES:
   use ncdf_river, only: init_river_input_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
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
   write(debug,*) 'init_river_input() # ',ncall
#endif

   LEVEL2 'init_river_input'

   select case (NETCDF)
      case (ANALYTICAL)
         LEVEL3 'Analytical boundary formulations'
      case (ASCII)
         LEVEL3 'ASCII boundary format'
      case (NETCDF)
         call init_river_input_ncdf(fn,n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'init_river_input'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_river_input()'
   write(debug,*)
#endif
   return
   end subroutine init_river_input
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
