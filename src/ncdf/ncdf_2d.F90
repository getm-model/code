!$Id: ncdf_2d.F90,v 1.2 2003-04-23 11:53:24 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_2d() - saves 2D-fields.
!
! !INTERFACE:
   module ncdf_2d
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

   integer                             :: x_dim,y_dim,time_dim

   integer                             :: ioff_id,joff_id,grid_type_id
   integer                             :: xc_id,xx_id,xu_id,xv_id
   integer                             :: yc_id,yx_id,yu_id,yv_id
   integer                             :: dx_id,dy_id
   integer                             :: lonc_id,latc_id
   integer                             :: time_id
   integer                             :: bathymetry_id
   integer                             :: elev_id,u_id,v_id
   integer                             :: res_u_id,res_v_id,surfdiv_id
   integer                             :: u10_id,v10_id
   integer                             :: airp_id,t2_id,hum_id,cc_id
   integer                             :: tausx_id,tausy_id,swr_id,shf_id
   integer                             :: xlen,ylen
   integer, parameter                  :: size_2d=150000
   REAL_4B                             :: ws(size_2d)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d.F90,v $
!  Revision 1.2  2003-04-23 11:53:24  kbk
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
