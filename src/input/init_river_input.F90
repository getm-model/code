!$Id: init_river_input.F90,v 1.2 2003-04-07 12:58:21 kbk Exp $
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
   character(len=*), intent(in) :: fn
   integer, intent(in)		:: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: init_river_input.F90,v $
!  Revision 1.2  2003-04-07 12:58:21  kbk
!  parallel + cleaned code
!
!  Revision 1.1.1.1  2002/05/02 14:01:34  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/10/07 14:50:22  bbh
!  Reading river data implemented - NetCFD
!
! !LOCAL VARIABLES:
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
