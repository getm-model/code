!$Id: ncdf_2d.F90,v 1.7 2008-09-16 11:21:50 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate 2D netCDF quantities
!
! !INTERFACE:
   module ncdf_2d
!
! !DESCRIPTION:
!
! !USES:
   use output
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: ncid=-1

   integer                             :: x_dim,y_dim
   integer                             :: time_dim
   integer                             :: time_id

   integer                             :: elev_id,u_id,v_id
   integer                             :: res_u_id,res_v_id,surfdiv_id
   integer                             :: u10_id,v10_id
   integer                             :: airp_id,t2_id,hum_id,tcc_id
   integer                             :: tausx_id,tausy_id,swr_id,shf_id
   integer                             :: evap_id=-1,precip_id=-1
   integer                             :: break_stat_id=-1

   REAL_4B, dimension(:), allocatable :: ws

! !DEFINED PARAMETERS
   REALTYPE, parameter                 :: elev_missing       =-9999.0
   REALTYPE, parameter                 :: vel_missing        =-9999.0
   REALTYPE, parameter                 :: airp_missing       =-9999.0
   REALTYPE, parameter                 :: t2_missing         =-9999.0
   REALTYPE, parameter                 :: hum_missing        =-9999.0
   REALTYPE, parameter                 :: tcc_missing        =-9999.0
   REALTYPE, parameter                 :: stress_missing     =-9999.0
   REALTYPE, parameter                 :: swr_missing        =-9999.0
   REALTYPE, parameter                 :: shf_missing        =-9999.0
   REALTYPE, parameter                 :: divergence_missing =-9999.0
   REALTYPE, parameter                 :: evap_missing       =-9999.0
   REALTYPE, parameter                 :: precip_missing     =-9999.0 
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d.F90,v $
!  Revision 1.7  2008-09-16 11:21:50  kb
!  if -DUSE_BREAKS save break statistics
!
!  Revision 1.6  2007-06-27 08:39:37  kbk
!  support for fresh water fluxes at the sea surface - Adolf Stips
!
!  Revision 1.5  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.4  2003/06/17 14:53:29  kbk
!  default meteo variables names comply with Adolf Stips suggestion + southpole(3)
!
!  Revision 1.3  2003/05/09 11:38:26  kbk
!  added proper undef support - based on Adolf Stips patch
!
!  Revision 1.2  2003/04/23 11:53:24  kbk
!  save lat/lon info for spherical grid
!
!  Revision 1.1.1.1  2002/05/02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/10/26 12:18:06  bbh
!  No actual storing of data in init_2d_ncdf.F90 -> save_2d_ncdf.F90
!
!  Revision 1.2  2001/09/27 08:35:10  bbh
!  Saving meteo again - in .2d.nc file
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
