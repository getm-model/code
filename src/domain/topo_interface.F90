!$Id: topo_interface.F90,v 1.1 2005-04-25 09:32:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: topo_interface --- for type safety of grid-related functions calls
!
! !INTERFACE:
   module topo_interface
!
! !DESCRIPTION:
! This module contains explicit interfaces for some external functions related 
! to checking and reading the topography and grid. These explicit interfaces
! prevent functions in the module {\tt domain} from calling the routines with
! any incorrect type and number of arguments (type-safety).
! 
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: topo_interface.F90,v $
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine check_grid(filename,filetype,iextr,jextr)
         character(len=*)              :: filename
         integer, intent(in)           :: filetype
#ifdef STATIC 
         integer, intent(in)           :: iextr
         integer, intent(in)           :: jextr
#else
         integer, intent(out)          :: iextr
         integer, intent(out)          :: jextr
#endif
      end subroutine check_grid
   end interface

   interface
      subroutine get_grid(filetype,H,Hland, &
                          iextr,jextr,ioff,joff,imin,imax,jmin,jmax)
         integer, intent(in)           :: filetype
         integer, intent(in)           :: iextr,jextr,ioff,joff
         integer, intent(in)           :: imin,imax,jmin,jmax

         REALTYPE, intent(out)         :: H(E2DFIELD)
         REALTYPE, intent(out)         :: Hland
      end subroutine get_grid
   end interface

!-----------------------------------------------------------------------

   end module topo_interface

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------


