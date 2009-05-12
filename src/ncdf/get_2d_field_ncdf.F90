!$Id: get_2d_field_ncdf.F90,v 1.2 2009-05-12 10:50:44 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: field_2d_ncdf - Interface with 2D field from file
!
! !INTERFACE:
   module field_2d_ncdf

!
! !DESCRIPTION:
! This module is responsible for reading 2D field quantities
! contained in a netCDF file. 
! !USES:
   use netcdf
   use exceptions
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
public get_2d_field_ncdf
!
! !PUBLIC DATA MEMBERS:
!
! !DEFINED PARAMETERS:
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
!

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: get_2d_field_ncdf()
!
! !INTERFACE:
   subroutine get_2d_field_ncdf(fn,varname,il,ih,jl,jh,field)
! !USES:
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
   character(len=*), intent(in)        :: fn,varname
   integer,          intent(in)        :: il,ih,jl,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(:,:)
!
! !LOCAL VARIABLES:
   integer, dimension(2)               :: start
   integer, dimension(2)               :: edges
   integer, dimension(2)               :: ubounds
   integer                             :: status,ncid,varid
!EOP
!-------------------------------------------------------------------------
!#include"netcdf.inc"

   LEVEL3 'get_2d_field_ncdf()'

   start(1) = il
   start(2) = jl
   edges(1) = ih-il+1
   edges(2) = jh-jl+1

   ubounds =  ubound(field)

   if ((ubounds(1) .ne. edges(1)) .or. ubounds(2) .ne. edges(2) ) then
      call getm_error("get_2d_field_ncdf()", &
           "Array bounds inconsistent.")
   endif

   status = nf90_open(trim(fn),NF90_NOWRITE,ncid)
   if (status .NE. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error opening file "//trim(fn))
   end if

   status = nf90_inq_varid(ncid,trim(varname),varid)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error inquiring "//trim(varname))
   endif

   status = nf90_get_var(ncid,varid,field,start,edges)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error reading "//trim(varname))
   endif

  status = nf90_close(ncid)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error closing file")
   endif

   return
   end subroutine get_2d_field_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Hans Burchard and Karsten Bolding (BB)          !
!-----------------------------------------------------------------------

end module field_2d_ncdf
