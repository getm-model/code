#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: get_2d_field_ncdf()
!
! !INTERFACE:
   subroutine get_2d_field_ncdf_by_id(ncid,varid,il,ih,jl,jh,n,field)
! !USES:
   use netcdf
   use exceptions
   IMPLICIT NONE
!
! !DESCRIPTION:
!  A two-dimensional netCDF variable with specified global range
! {\tt il < i < ih} and {\tt jl < j < jh} is read into {\tt field}.
! It is checked if the sizes of the fields correspond exactly.
! When calling this funtions, remember that  FORTRAN netCDF variables
! start with index 1.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,varid
   integer, intent(in)                 :: il,ih,jl,jh,n
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Lars Umlauf
!
! !LOCAL VARIABLES:
   integer, dimension(3)               :: start
   integer, dimension(3)               :: edges
   integer, dimension(2)               :: ubounds
   integer                             :: status
!EOP
!-------------------------------------------------------------------------
!KB   LEVEL3 'get_2d_field_ncdf_by_id()'
   start(1) = il
   start(2) = jl
   start(3) = n
   edges(1) = ih-il+1
   edges(2) = jh-jl+1
   edges(3) = 1

   ubounds =  ubound(field)

   if ((ubounds(1) .ne. edges(1)) .or. ubounds(2) .ne. edges(2) ) then
      call getm_error("get_2d_field_ncdf_by_id()", &
           "Array bounds inconsistent.")
   endif

   status = nf90_get_var(ncid,varid,field,start,edges)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf_by_id()", &
           "Error reading variable")
   endif

   return
   end subroutine get_2d_field_ncdf_by_id
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding (BB)          !
!-----------------------------------------------------------------------
