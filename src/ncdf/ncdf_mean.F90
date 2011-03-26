#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate netCDF mean quantities
!
! !INTERFACE:
   module ncdf_mean
!
! !DESCRIPTION:
!
! !USES:
   use output
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: ncid=-1

   integer                             :: x_dim,y_dim,z_dim
   integer                             :: time_dim
   integer                             :: time_id

   integer                             :: swrmean_id,ustarmean_id,ustar2mean_id
   integer                             :: uumean_id,vvmean_id,wmean_id
   integer                             :: saltmean_id,tempmean_id,hmean_id
   integer                             :: nm3dS_id,nm3dT_id,nm2dS_id,nm2dT_id
   integer                             :: pm3dS_id,pm3dT_id,pm2dS_id,pm2dT_id
#ifdef GETM_BIO
   integer, allocatable                :: biomean_id(:)
#endif

   REALTYPE, parameter                 :: hh_missing=-10.0
   REALTYPE, parameter                 :: swr_missing=-9999.0
   REALTYPE, parameter                 :: vel_missing=-9999.0
   REALTYPE, parameter                 :: salt_missing=-9999.0
   REALTYPE, parameter                 :: temp_missing=-9999.0
   REALTYPE, parameter                 :: tke_missing=-9999.0
   REALTYPE, parameter                 :: eps_missing=-9999.0
   REALTYPE, parameter                 :: nummix_missing=-9999.0
#ifdef GETM_BIO
   REALTYPE, parameter                 :: bio_missing=-9999.0
#endif

!
!  Original author(s): Adolf Stips & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_mean

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
