#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: read bathymetry and grid information
!
! !INTERFACE:
   subroutine read_topo_file(filetype,grid_type)
!
! !DESCRIPTION:
!  This is a wrapper routine to read the bathymetry file
!  {\tt filename} containing the bathymetry and grid infformation.
!  The only thing it actually does is calling the specialised routines
!  according to the {\tt filetype}. Presently only NetCDF format is
!  supported
!
! !USES:
  use exceptions
  use ncdf_topo, only : ncdf_read_topo_file
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: filetype
   integer,intent(in)                  :: grid_type
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 "Reading bathymetry and grid information ..."

   select case (filetype)
      case(ASCII)
         call getm_error("get_grid()", &
              "ASCII format for topo files not yet supported.")
      case(NETCDF)
         call ncdf_read_topo_file(grid_type)
      case(BINARY)
         call getm_error("get_grid()", &
              "RAWBINDARY format for topo files not yet supported.")
      case default
         call getm_error("get_grid()","unkown format for topo file.")
   end select

   LEVEL2 "done - reading bathymetry and grid"

   return
   end subroutine read_topo_file
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------

