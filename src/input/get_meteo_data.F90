!$Id: get_meteo_data.F90,v 1.2 2003-04-23 12:04:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_meteo_data - read meteo data from file.
!
! !INTERFACE:
   subroutine get_meteo_data(n)
!
! !DESCRIPTION:
!  Reads specified boundary data from file(s)
!
! !USES:
   use ncdf_meteo, only: get_meteo_data_ncdf
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
!  $Log: get_meteo_data.F90,v $
!  Revision 1.2  2003-04-23 12:04:08  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:33  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.1  2001/05/25 19:12:51  bbh
!  ncdf/ncdf_meteo.F90
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'get_meteo_data() # ',ncall
#endif

   select case (NETCDF)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'should get ASCII boundary data'
      case (NETCDF)
         call get_meteo_data_ncdf(n)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_meteo_data'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_meteo_data()'
   write(debug,*)
#endif
   return
   end subroutine get_meteo_data
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
