!$Id: ncdf_3d.F90,v 1.9 2009-01-05 09:57:06 kb Exp $
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
   integer                             :: uu_id,vv_id,w_id
   integer                             :: salt_id,temp_id,sigma_t_id
   integer                             :: rad_id
   integer                             :: tke_id,num_id,nuh_id,eps_id
   integer                             :: SS_id,NN_id
#ifdef SPM
   integer                             :: spmpool_id,spm_id
#endif
#ifdef GETM_BIO
   integer, allocatable                :: bio_ids(:)
#endif

   REAL_4B, dimension(:), allocatable  :: ws

! !DEFINED PARAMETERS
   REALTYPE, parameter                 :: hh_missing     =-9999.0
   REALTYPE, parameter                 :: elev_missing   =-9999.0
   REALTYPE, parameter                 :: vel_missing    =-9999.0
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
#ifdef SPM
   REALTYPE, parameter                 :: spmpool_missing=-9999.0
   REALTYPE, parameter                 :: spm_missing    =-9999.0
#endif
#ifdef GETM_BIO
   REALTYPE, parameter                 :: bio_missing=-9999.0
#endif

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d.F90,v $
!  Revision 1.9  2009-01-05 09:57:06  kb
!  option for storing SS and NN
!
!  Revision 1.8  2007-02-20 13:52:15  kbk
!  solar radiation -> 3d field - possible to save
!
!  Revision 1.7  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.6  2004/06/15 08:25:57  kbk
!  added supoort for spm - Ruiz
!
!  Revision 1.5  2004/05/04 09:23:51  kbk
!  hydrostatic consistency criteria stored in .3d.nc file
!
!  Revision 1.4  2003/12/16 12:51:04  kbk
!  preparing for proper support for SPM (manuel)
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
!  Revision 1.3  2001/10/23 14:19:20  bbh
!  Stores h if general vertical coordinates
!
!  Revision 1.2  2001/10/23 07:37:17  bbh
!  Saving spm - if calc_spm and save_spm are both true
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
