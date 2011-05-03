#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate 3D netCDF quantities
!
! !INTERFACE:
   module ncdf_3d
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

   integer                             :: hcc_id,h_id
   integer                             :: elev_id,u_id,v_id
   integer                             :: taubx_id,tauby_id
   integer                             :: uu_id,vv_id,w_id
   integer                             :: salt_id,temp_id,sigma_t_id
   integer                             :: rad_id
   integer                             :: tke_id,num_id,nuh_id,eps_id
   integer                             :: SS_id,NN_id
#ifdef _LES_
   integer                             :: Am3d_id
#endif
#ifdef SPM
   integer                             :: spmpool_id,spm_id
#endif
#ifdef GETM_BIO
   integer, allocatable                :: bio_ids(:)
#endif
#ifdef _FABM_
   integer, allocatable, dimension(:)  :: fabm_ids,fabm_ids_diag,fabm_ids_ben,fabm_ids_diag_hz
#endif
   integer                             :: nm3dS_id,nm3dT_id,nm2dS_id,nm2dT_id
   integer                             :: pm3dS_id,pm3dT_id,pm2dS_id,pm2dT_id

   REALTYPE, dimension(:,:,:), allocatable :: ws

! !DEFINED PARAMETERS
   REALTYPE, parameter                 :: hh_missing     =-9999.0
   REALTYPE, parameter                 :: elev_missing   =-9999.0
   REALTYPE, parameter                 :: vel_missing    =-9999.0
   REALTYPE, parameter                 :: tau_missing    =-9999.0
   REALTYPE, parameter                 :: salt_missing   =-9999.0
   REALTYPE, parameter                 :: temp_missing   =-9999.0
   REALTYPE, parameter                 :: rho_missing    =-9999.0
   REALTYPE, parameter                 :: rad_missing    =-9999.0
   REALTYPE, parameter                 :: tke_missing    =-9999.0
   REALTYPE, parameter                 :: nuh_missing    =-9999.0
   REALTYPE, parameter                 :: num_missing    =-9999.0
   REALTYPE, parameter                 :: eps_missing    =-9999.0
   REALTYPE, parameter                 :: SS_missing     =-9999.0
   REALTYPE, parameter                 :: NN_missing     =-9999.0
#ifdef _LES_
   REALTYPE, parameter                 :: Am3d_missing   =-9999.0
#endif
#ifdef SPM
   REALTYPE, parameter                 :: spmpool_missing=-9999.0
   REALTYPE, parameter                 :: spm_missing    =-9999.0
#endif
#if (defined(GETM_BIO) || defined(_FABM_))
   REALTYPE, parameter                 :: bio_missing=-9999.0
#endif
   REALTYPE, parameter                 :: nummix_missing=-9999.0

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
