!$Id: get_river_data.F90,v 1.4 2003-12-16 12:33:38 kbk Exp $
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: get_river_data.F90,v $
!  Revision 1.4  2003-12-16 12:33:38  kbk
!  some typos (manuel)
!
!  Revision 1.3  2003/04/23 12:04:08  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:58:21  kbk
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
