!$Id: ncdf_mean.F90,v 1.1 2004-03-29 15:38:10 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_mean() - saves Mean-fields.
!
! !INTERFACE:
   module ncdf_mean
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
   integer         :: ncid=-1

   integer         :: x_dim,y_dim,z_dim,time_dim

   integer         :: ioff_id,joff_id,grid_type_id,vert_cord_id
   integer         :: xc_id,xx_id,xu_id,xv_id
   integer         :: yc_id,yx_id,yu_id,yv_id
   integer         :: z_id
   integer         :: dx_id,dy_id
   integer         :: lonc_id,latc_id
   integer         :: time_id
   integer         :: bathymetry_id
   integer         :: h_id=-1
   integer         :: swrmean_id,ustarmean_id,ustar2mean_id
   integer         :: uumean_id,vvmean_id,wmean_id
   integer         :: saltmean_id,tempmean_id,hmean_id
   integer         :: xlen,ylen,zlen
   integer, parameter        :: size_3d=9000000
   REAL_4B         :: ws(size_3d)
   REALTYPE, parameter       :: h_missing=-10.0
   REALTYPE, parameter       :: swr_missing=-9999.0
   REALTYPE, parameter       :: vel_missing=-9999.0
   REALTYPE, parameter       :: salt_missing=-9999.0
   REALTYPE, parameter       :: temp_missing=-9999.0
   REALTYPE, parameter       :: tke_missing=-9999.0
   REALTYPE, parameter       :: eps_missing=-9999.0
!
!  Original author(s): Adolf Stips & Karsten Bolding
!
!  $Log: ncdf_mean.F90,v $
!  Revision 1.1  2004-03-29 15:38:10  kbk
!  possible to store calculated mean fields
!
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_mean

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
