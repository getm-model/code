#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_2d_bdy - read boundary data file(s)
!
! !INTERFACE:
   subroutine get_2d_bdy(fmt,n)
!
! !DESCRIPTION:
!  Reads specified boundary data from file(s)
!
! !USES:
   use ncdf_2d_bdy, only: do_2d_bdy_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: fmt,n
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
   write(debug,*) 'get_2d_bdy() # ',ncall
#endif
   select case (fmt)
      case (NO_DATA)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'should get ASCII boundary data'
      case (NETCDF)
         call do_2d_bdy_ncdf(n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_2d_bdy'
      end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_2d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine get_2d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
