#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_3d_bdy - initialise 3D boundary data file(s)
!
! !INTERFACE:
   subroutine init_3d_bdy(fn,fmt,n)
!
! !DESCRIPTION:
!  Prepares for reading boundary data for the 2D module.
!
! !USES:
   use ncdf_3d_bdy, only: init_3d_bdy_ncdf
#ifdef _FABM_
   use gotm_fabm, only: fabm_calc
   use ncdf_3d_bio_bdy, only: init_3d_bio_bdy_ncdf
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: fmt,n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: bdyfmt=NETCDF
#ifdef _FABM_
   character(len=255)        :: bio_fn
#endif
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
         call init_3d_bdy_ncdf(fn,n)
#ifdef _FABM_
         if (fabm_calc) then
            bio_fn='bdy_3d_bio.nc'
            LEVEL3 'reading BIO from: ',trim(bio_fn)
            call init_3d_bio_bdy_ncdf(bio_fn)
         end if
#endif
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
