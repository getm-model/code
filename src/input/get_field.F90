#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: get_field()
!
! !INTERFACE:
   module get_field
!
! !DESCRIPTION:
!  Provides generic subroutines for reading 2D and 3D fields. Presently
!  only NetCDF is supported.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public get_2d_field, get_3d_field
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_2d_field - read a 2D field from a file.
!
! !INTERFACE:
   subroutine get_2d_field(fn,varname,il,ih,jl,jh,break_on_missing,field)
!
! !DESCRIPTION:
!  Reads varname from a named file - fn - into to field.
!
! !USES:
   use ncdf_get_field, only: get_2d_field_ncdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,varname
   integer, intent(in)                 :: il,ih,jl,jh
   logical, intent(in)                 :: break_on_missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer, parameter        :: fmt=NETCDF
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   ncall = ncall+1
   write(debug,*) 'get_field() # ',ncall
#endif

   select case (fmt)
      case (ANALYTICAL)
      case (ASCII)
         STDERR 'Should get an ASCII field'
         stop 'get_2d_field()'
      case (NETCDF)
         call get_2d_field_ncdf(fn,varname,il,ih,jl,jh,break_on_missing,field)
      case DEFAULT
         FATAL 'A non valid input format has been chosen'
         stop 'get_2d_field'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving get_2d_field()'
   write(debug,*)
#endif
   return
   end subroutine get_2d_field
!EOC

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

   end module get_field

!-----------------------------------------------------------------------
! Copyright (C) 2015 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
