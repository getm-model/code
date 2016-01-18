#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ice_input - initialise 2D boundary data file(s)
!
! !INTERFACE:
   subroutine init_ice_input(fn,n)
!
! !DESCRIPTION:
!  Prepares for reading boundary data for the 2D module.
!
! !USES:
   use ncdf_ice, only: init_ice_input_ncdf
!kbk   use ice, only: airp,tausx,tausy
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
   write(debug,*) 'init_ice_input() # ',ncall
#endif

   LEVEL2 'init_ice_input'

   select case (NETCDF)
      case (ANALYTICAL)
         LEVEL3 'Analytical boundary formulations'
      case (ASCII)
         LEVEL3 'ASCII boundary format'
      case (NETCDF)
         call init_ice_input_ncdf(fn,n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'init_ice_input'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_ice_input()'
   write(debug,*)
#endif
   return
   end subroutine init_ice_input
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
