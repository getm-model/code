!$Id: init_2d_bdy.F90,v 1.1 2002-05-02 14:01:33 gotm Exp $
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
   character(len=*), intent(in)	:: fn
   integer, intent(in)	:: fmt
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: init_2d_bdy.F90,v $
!  Revision 1.1  2002-05-02 14:01:33  gotm
!  Initial revision
!
!  Revision 1.5  2001/09/19 11:20:32  bbh
!  Explicit de-allocates memory when -DFORTRAN90
!
!  Revision 1.4  2001/05/25 19:31:45  bbh
!  Removed LEVEL3 statement
!
!  Revision 1.3  2001/05/14 12:45:24  bbh
!  Introduced module ncdf_2d_bdy
!
!  Revision 1.2  2001/05/10 11:35:41  bbh
!  Now calls init_2d_bdy_ncdf()
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
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
