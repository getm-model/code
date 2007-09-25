!$Id: get_grid.F90,v 1.3 2007-09-25 08:31:18 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: read bathymetry
!
! !INTERFACE:
   subroutine get_grid(filetype,H,Hland,iextr,jextr,ioff,joff,imin,imax,jmin,jmax)
!
! !DESCRIPTION:
!  This is a wrapper routine to read the bathymetry file
!  {\tt filename} containing the grid and the bathymetery. {\tt get\_grid}
!  assumes that {\tt check\_grid} has been called before.
!  The only thing it actually does is calling the specialised routines 
!  according to the {\tt filetype}.
!
! !USES:
  use exceptions
  use ncdf_topo, only : ncdf_get_grid
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: filetype
   integer, intent(in)                 :: iextr,jextr,ioff,joff
   integer, intent(in)                 :: imin,imax,jmin,jmax
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: H(E2DFIELD)
   REALTYPE, intent(out)               :: Hland
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: get_grid.F90,v $
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

   LEVEL2 "Reading grid and bathymetry..."

   select case (filetype)
      case(ASCII)
         call getm_error("get_grid()", &
              "ASCII format for topo files not yet supported.")
      case(NETCDF)
         call ncdf_get_grid(H,Hland,iextr,jextr,ioff,joff,imin,imax,jmin,jmax)
      case(BINARY)
         call getm_error("get_grid()", &
              "RAWBINDARY format for topo files not yet supported.")
      case default
         call getm_error("get_grid()","unkown format for topo file.")
   end select

   LEVEL2 "...done"

   return
   end subroutine get_grid
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------

