#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_3d_field - read a 3D field from a file.
!
! !INTERFACE:
   subroutine get_3d_field(fname,var,n,break_on_missing,f)
!
! !DESCRIPTION:
!  Reads from file - fname - the variable var into f.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use ncdf_get_field, only: get_3d_field_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname,var
   integer, intent(in)                 :: n
   logical, intent(in)                 :: break_on_missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer, parameter        :: fmt=NETCDF
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'get_3d_field() # ',ncall
#endif
   select case (fmt)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'Should get an ASCII field'
         stop 'get_3d_field()'
      case (NETCDF)
         call get_3d_field_ncdf(fname,var,n,break_on_missing,f)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_3d_field'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_3d_field()'
   write(debug,*)
#endif
   return
   end subroutine get_3d_field
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
