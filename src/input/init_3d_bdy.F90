!$Id: init_3d_bdy.F90,v 1.1 2002-05-02 14:01:35 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_3d_bdy - initialise 3D boundary data file(s)
!
! !INTERFACE:
   subroutine init_3d_bdy(fn,fmt)
!
! !DESCRIPTION:
!  Prepares for reading boundary data for the 2D module.
!
! !USES:
   use ncdf_3d_bdy, only: init_3d_bdy_ncdf
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
!  $Log: init_3d_bdy.F90,v $
!  Revision 1.1  2002-05-02 14:01:35  gotm
!  Initial revision
!
!
! !LOCAL VARIABLES:
   integer	:: rc
   integer	:: bdyfmt=NETCDF
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'init_3d_bdy() # ',ncall
#endif

   LEVEL2 'init_3d_bdy'

   select case (fmt)
      case (ANALYTICAL)
         LEVEL3 'Analytical boundary formulations'
         stop 'init_3d_bdy'
      case (ASCII)
         LEVEL3 'ASCII boundary format'
         stop 'init_3d_bdy'
      case (NETCDF)
         LEVEL3 'reading from: ',trim(fn)
         call init_3d_bdy_ncdf(fn)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'init_3d_bdy'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_3d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine init_3d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
