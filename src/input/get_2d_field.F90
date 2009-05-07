!$Id: get_2d_field.F90,v 1.1 2009-05-07 16:02:01 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_2d_field - read a 2D field from a file.
!
! !INTERFACE:
   subroutine get_2d_field(fn,varname,il,ih,jl,jh,field)
!
! !DESCRIPTION:
!  Reads varname from a named file - fn - into to field.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,varname
   integer, intent(in)                 :: il,ih,jl,jh
!
! !INPUT/OUTPUT PARAMETERS:
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
         call get_2d_field_ncdf(fn,varname,il,ih,jl,jh,field)
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
! Copyright (C) 2009 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
