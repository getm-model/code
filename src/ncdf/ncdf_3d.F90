!$Id: ncdf_3d.F90,v 1.3 2003-05-09 11:38:26 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_3d() - saves 2D-fields.
!
! !INTERFACE:
   module ncdf_3d
!
! !DESCRIPTION:
!
! !USES:
   use output
   use ncdf_common
   IMPLICIT NONE
!
!   private
!
! !PUBLIC DATA MEMBERS:
   integer, public                     :: ncid=-1

   integer                             :: x_dim,y_dim,z_dim,time_dim

   integer                             :: ioff_id,joff_id,grid_type_id,vert_cord_id
   integer                             :: xc_id,xx_id,xu_id,xv_id
   integer                             :: yc_id,yx_id,yu_id,yv_id
   integer                             :: z_id
   integer                             :: dx_id,dy_id
   integer                             :: lonc_id,latc_id
   integer                             :: time_id
   integer                             :: bathymetry_id
   integer                             :: h_id=-1
   integer                             :: elev_id,u_id,v_id
   integer                             :: uu_id,vv_id,w_id
   integer                             :: salt_id,temp_id,sigma_t_id
   integer                             :: tke_id,num_id,nuh_id,eps_id
   integer                             :: spm_id

   integer                             :: xlen,ylen,zlen
   integer, parameter                  :: size_3d=9000000
   REAL_4B                             :: ws(size_3d)
   REALTYPE, parameter                 :: h_missing=-10.0
   REALTYPE, parameter                 :: elev_missing=-9999.0
   REALTYPE, parameter                 :: vel_missing=-9999.0
   REALTYPE, parameter                 :: salt_missing=-9999.0
   REALTYPE, parameter                 :: temp_missing=-9999.0
   REALTYPE, parameter                 :: rho_missing=-9999.0
   REALTYPE, parameter                 :: tke_missing=-9999.0
   REALTYPE, parameter                 :: nuh_missing=-9999.0
   REALTYPE, parameter                 :: num_missing=-9999.0
   REALTYPE, parameter                 :: eps_missing=-9999.0
   REALTYPE, parameter                 :: spm_missing=-9999.0


!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d.F90,v $
!  Revision 1.3  2003-05-09 11:38:26  kbk
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
