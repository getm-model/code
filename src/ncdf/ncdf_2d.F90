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
   integer                             :: res_u_id,res_v_id
   integer                             :: u10_id,v10_id
   integer                             :: airp_id,t2_id,hum_id,tcc_id
   integer                             :: tausx_id,tausy_id,swr_id,shf_id
   integer                             :: evap_id=-1,precip_id=-1
   integer                             :: Am_2d_id
   integer                             :: break_stat_id=-1

   REALTYPE, dimension(:,:), allocatable :: ws

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
   REALTYPE, parameter                 :: evap_missing       =-9999.0
   REALTYPE, parameter                 :: precip_missing     =-9999.0
   REALTYPE, parameter                 :: Am_2d_missing      =-9999.0
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
