!$Id: check_grid.F90,v 1.3 2007-09-25 08:31:18 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: check_grid - generic subroutine for reading the grid
!
! !INTERFACE:
   subroutine check_grid(filename,filetype,iextr,jextr)
!
! !DESCRIPTION:
!  This is a wrapper routine to check whether the bathymetry file
!  {\tt filename} contains everything needed for the desired 
!  {\tt grid\_type}. For {\tt STATIC} allocation, it checks whether 
!  {\tt iextr} and {\tt jextr} are in agreement with the size of the 
!  read grid. For dynamic  allocation, it derives {\tt iextr} and 
!  {\tt jextr} from the dimensions of the grid.
!  The only thing {\tt check\_grid} actually does is calling specialised 
!  routines according to the {\tt filetype} of the file {\tt filename}.
!  
! !USES:
  use exceptions
  use ncdf_topo, only : ncdf_check_grid
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)                :: filename
   integer, intent(in)             :: filetype
#ifdef STATIC 
   integer, intent(in)             :: iextr
   integer, intent(in)             :: jextr
#else
   integer, intent(out)            :: iextr
   integer, intent(out)            :: jextr
#endif
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: check_grid.F90,v $
!  Revision 1.3  2007-09-25 08:31:18  kbk
!  RAWBINARY --> BINARY
!
!  Revision 1.2  2006-01-29 20:32:34  hb
!  Small LaTeX corrections to source code documentation
!
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 "Checking grid and bathymetry..."

   select case (filetype)
      case(ASCII)
         call getm_error("check_grid()", &
              "ASCII format for "//trim(filename)//" not yet supported.")
      case(NETCDF)
         call ncdf_check_grid(filename,iextr,jextr)
      case(BINARY)
         call getm_error("check_grid()", &
              "RAWBINDARY format for "//trim(filename)//" not yet supported.")
      case default
         call getm_error("check_grid()", &
              "unkown format for"//trim(filename)//".")
   end select

   LEVEL2 "...done"

   return
   end subroutine check_grid
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------

