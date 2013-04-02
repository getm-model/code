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

   integer                             :: hcc_id,h_id=-1
   integer                             :: elev_id
   integer                             :: fluxu_adv_id=-1
   integer                             :: fluxv_adv_id=-1
   integer                             :: u_adv_id=-1
   integer                             :: v_adv_id=-1
   integer                             :: fluxuu_id=-1
   integer                             :: fluxvv_id=-1
   integer                             :: fluxw_id=-1
   integer                             :: uu_id=-1
   integer                             :: vv_id=-1
   integer                             :: w_id=-1
   integer                             :: taubx_id,tauby_id

#ifdef _MOMENTUM_TERMS_
   integer                             :: tdv_u_id
   integer                             :: adv_u_id
   integer                             :: vsd_u_id
   integer                             :: hsd_u_id
   integer                             :: cor_u_id
   integer                             :: epg_u_id
   integer                             :: ipg_u_id

   integer                             :: tdv_v_id
   integer                             :: adv_v_id
   integer                             :: vsd_v_id
   integer                             :: hsd_v_id
   integer                             :: cor_v_id
   integer                             :: epg_v_id
   integer                             :: ipg_v_id
#endif
#if defined(CURVILINEAR)
   integer                             :: uurot_id,vvrot_id
#endif
   integer                             :: salt_id,temp_id,sigma_t_id
   integer                             :: rad_id
   integer                             :: tke_id,num_id,nuh_id,eps_id
   integer                             :: SS_id,NN_id
   integer                             :: Am_3d_id=-1
   integer                             :: diffxx_id=-1
   integer                             :: diffyy_id=-1
   integer                             :: diffxy_id=-1
#ifdef SPM
   integer                             :: spmpool_id,spm_id
#endif
#ifdef GETM_BIO
   integer, allocatable                :: bio_ids(:)
#endif
#ifdef _FABM_
   integer, allocatable, dimension(:)  :: fabm_ids,fabm_ids_diag,fabm_ids_ben,fabm_ids_diag_hz
   integer, allocatable, dimension(:)  :: nmpel_ids
#endif
   integer                             :: nd3d_id=-1
   integer                             :: nd3do_id=-1,pd3d_id=-1
   integer                             :: nmS_id=-1
   integer                             :: nmT_id=-1
   integer                             :: nmSo_id=-1,pmS_id=-1
   integer                             :: nmTo_id=-1,pmT_id=-1

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
   REALTYPE, parameter                 :: Am_3d_missing  =-9999.0
   REALTYPE, parameter                 :: stirr_missing  =-9999.0
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
