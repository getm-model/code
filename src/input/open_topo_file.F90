#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: open_topo_file - generic subroutine for reading the grid
!
! !INTERFACE:
   subroutine open_topo_file(filetype,filename,grid_type,iextr,jextr)
!
! !DESCRIPTION:
!  This is a wrapper routine for calling specialised routines according
!  to the {\tt filetype} of the bathymetry file {\tt filename} to open
!  it and to inquire {\tt grid\_type}, {\tt iextr} and {\tt jextr}.
!  For {\tt STATIC} allocation, the consistency check for
!  {\tt iextr} and {\tt jextr} has to be done in the calling routine!
!
! !USES:
  use exceptions
  use ncdf_topo, only : ncdf_open_topo_file
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)             :: filetype
   character(len=*)                :: filename
!
! !OUTPUT PARAMETERS:
   integer, intent(out)            :: grid_type
   integer, intent(out)            :: iextr,jextr
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
!BOC

   select case (filetype)
      case(ASCII)
         call getm_error("open_topo_file()", &
              "ASCII format for "//trim(filename)//" not yet supported.")
      case(NETCDF)
         call ncdf_open_topo_file(filename,grid_type,iextr,jextr)
      case(BINARY)
         call getm_error("open_topo_file()", &
              "RAWBINDARY format for "//trim(filename)//" not yet supported.")
      case default
         call getm_error("open_topo_file()", &
              "unkown format for"//trim(filename)//".")
   end select

   return
   end subroutine open_topo_file
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2014 - Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------

