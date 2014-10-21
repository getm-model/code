#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_waves_input -
!
! !INTERFACE:
   subroutine init_waves_input(fn,n)
!
! !DESCRIPTION:
!
! !USES:
   use ncdf_waves, only: init_waves_input_ncdf
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
   write(debug,*) 'init_waves_input() # ',ncall
#endif

   LEVEL2 'init_waves_input'

   select case (NETCDF)
      case (ANALYTICAL)
         LEVEL3 'Analytical boundary formulations'
      case (ASCII)
         LEVEL3 'ASCII boundary format'
      case (NETCDF)
         call init_waves_input_ncdf(fn,n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'init_waves_input'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_waves_input()'
   write(debug,*)
#endif
   return
   end subroutine init_waves_input
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2014 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
