#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_3d_bdy - read boundary data file(s)
!
! !INTERFACE:
   subroutine get_3d_bdy(fmt,n)
!
! !DESCRIPTION:
!  Reads specified boundary data from file(s)
!
! !USES:
   use ncdf_3d_bdy, only: do_3d_bdy_ncdf
   use m3d, only: calc_salt,calc_temp
#ifdef _FABM_
   use getm_fabm, only: fabm_calc
   use ncdf_3d_bio_bdy, only: do_3d_bio_bdy_ncdf
#endif
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
   write(debug,*) 'get_3d_bdy() # ',ncall
#endif
   select case (fmt)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'should get ASCII boundary data'
      case (NETCDF)
         if (calc_salt .or. calc_temp) then
            call do_3d_bdy_ncdf(n)
         end if
#ifdef _FABM_
         if (fabm_calc) then
            call do_3d_bio_bdy_ncdf(n)
         end if
#endif
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_3d_bdy'
      end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_3d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine get_3d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
