#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate grid related quantities
!
! !INTERFACE:
   module grid_ncdf
!
! !DESCRIPTION:
! This module is a container for grid related variables and
! parameters which are used jointly by different parts of the
! netCDF storage system.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: xlen=-1,ylen=-1,zlen=-1
   integer                             :: xc_dim=-1,yc_dim=-1
   integer                             :: xx_dim=-1,yx_dim=-1

! !DEFINED PARAMETERS
   REALTYPE, parameter                 :: h_missing      =-10.0
   REALTYPE, parameter                 :: xy_missing     =-999.0
   REALTYPE, parameter                 :: latlon_missing =-999.0
   REALTYPE, parameter                 :: conv_missing   =-999.0
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------

   end module grid_ncdf

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
