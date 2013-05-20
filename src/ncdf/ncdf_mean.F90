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
   integer                             :: hmean_id
   integer                             :: saltmean_id=-1
   integer                             :: tempmean_id=-1
   integer                             :: sigma_t_id=-1
   integer                             :: nm3dS_id,nm3dT_id,nm2dS_id,nm2dT_id
   integer                             :: pm3dS_id,pm3dT_id,pm2dS_id,pm2dT_id
   integer                             :: nm3d_id,nm2d_id
#ifdef GETM_BIO
   integer, allocatable                :: biomean_id(:)
#endif
#ifdef _FABM_
   integer, allocatable                :: fabmmean_ids(:)
   integer, allocatable                :: fabmmean_ids_ben(:)
   integer, allocatable                :: fabmmean_ids_diag(:)
   integer, allocatable                :: fabmmean_ids_diag_hz(:)
#endif

   REALTYPE, parameter                 :: hh_missing=-10.0
   REALTYPE, parameter                 :: swr_missing=-9999.0
   REALTYPE, parameter                 :: vel_missing=-9999.0
   REALTYPE, parameter                 :: salt_missing=-9999.0
   REALTYPE, parameter                 :: temp_missing=-9999.0
   REALTYPE, parameter                 :: rho_missing=-9999.0
   REALTYPE, parameter                 :: tke_missing=-9999.0
   REALTYPE, parameter                 :: eps_missing=-9999.0
   REALTYPE, parameter                 :: nummix_missing=-9999.0
#if (defined(GETM_BIO) || defined(_FABM_))
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
