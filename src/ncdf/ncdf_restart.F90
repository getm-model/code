#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate netCDF restart quantities
!
! !INTERFACE:
   module ncdf_restart
!
! !DESCRIPTION:
!  This module and the related *\_restart\_ncdf() subroutines provide a
!  drop-in replacement for the binary file hotstart facility in GETM.
!  The main reason for using NetCDF formatted hotstart files instead of
!  binary format is the abillity to use standard tools (nco, ncmerge) is
!  a much easier way to to introduce a new subdomain decomposition for
!  an already running set-up - without having to start all over again.
!  See $read\_restart\_ncdf()$ for further explanation.\newline
!  This modules just contains variables shared accros the *\_restart\_ncdf()
!  routines.
!
! !USES:
   use output
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: ncid=-1
   integer                             :: xdim_id=-1
   integer                             :: ydim_id=-1
   integer                             :: zdim_id=-1
   integer                             :: xax_id
   integer                             :: yax_id
   integer                             :: zax_id
   integer                             :: loop_id
   integer                             :: julianday_id
   integer                             :: secondsofday_id
   integer                             :: timestep_id
   integer                             :: z_id,zo_id
   integer                             :: U_id
   integer                             :: SlUx_id,Slru_id
   integer                             :: V_id
   integer                             :: SlVx_id,Slrv_id
#ifndef NO_3D
   integer                             :: ssen_id,ssun_id,ssvn_id
   integer                             :: sseo_id,ssuo_id,ssvo_id
   integer                             :: Uint_id,Vint_id
   integer                             :: uu_id,vv_id,ww_id
   integer                             :: uuEx_id,vvEx_id
   integer                             :: tke_id,eps_id
   integer                             :: num_id,nuh_id
   integer                             :: hn_id
#ifndef NO_BAROCLINIC
   integer                             :: T_id,S_id
#endif
   integer                             :: minus_bnh_id,wco_id
#ifdef SPM
   integer                             :: spm_id,spmpool_id
#endif
#ifdef GETM_BIO
   integer                             :: biodim_id
   integer                             :: bio_id
#endif
#ifdef _FABM_
   integer                             :: fabmpeldim_id
   integer                             :: fabmbendim_id
   integer                             :: fabm_pel_id
   integer                             :: fabm_ben_id
#endif
#endif

   integer                             :: xlen,ylen,zlen
   integer                             :: status
   integer                             :: start(5),edges(5)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_restart

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
