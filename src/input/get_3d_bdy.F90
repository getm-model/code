!$Id: get_3d_bdy.F90,v 1.1 2002-05-02 14:01:35 gotm Exp $
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: fmt,n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: get_3d_bdy.F90,v $
!  Revision 1.1  2002-05-02 14:01:35  gotm
!  Initial revision
!
!  Revision 1.3  2001/10/22 08:06:42  bbh
!  Removed a bogus subroutine call
!
!  Revision 1.2  2001/05/14 12:45:24  bbh
!  Introduced module ncdf_2d_bdy
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
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
         call do_3d_bdy_ncdf(n)
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
