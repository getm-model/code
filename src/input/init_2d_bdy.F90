#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_2d_bdy - initialise 2D boundary data file(s)
!
! !INTERFACE:
   subroutine init_2d_bdy(fn,fmt)
!
! !DESCRIPTION:
!  Prepares for reading boundary data for the 2D module.
!
! !USES:
   use ncdf_2d_bdy, only: init_2d_bdy_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: fmt
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
   write(debug,*) 'init_2d_bdy() # ',ncall
#endif

   LEVEL2 'init_2d_bdy'

   select case (fmt)
      case (NO_DATA)
      case (ANALYTICAL)
         LEVEL3 'Analytical boundary formulations'
         stop 'init_2d_bdy'
      case (ASCII)
         LEVEL3 'ASCII boundary format'
         stop 'init_2d_bdy'
      case (NETCDF)
         call init_2d_bdy_ncdf(fn)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'init_2d_bdy'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_2d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine init_2d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
