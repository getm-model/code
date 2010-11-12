!$Id: save_mean_ncdf.F90,v 1.3 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise mean netCDF variables
!
! !INTERFACE:
   subroutine save_mean_ncdf(secs)
!
! !DESCRIPTION:
!
! !USES:
   use exceptions
   use ncdf_mean
   use grid_ncdf
   use diagnostic_variables
   use domain,       only: ioff,joff,imin,imax,jmin,jmax,kmax
   use domain,       only: H,az
   use variables_3d, only: kmin
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: secs
!
! !DEFINED PARAMTERS:
   logical, parameter   :: save3d=.true.
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips & Karsten Bolding
!
!  $Log: save_mean_ncdf.F90,v $
!  Revision 1.3  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.2  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.1  2004/03/29 15:38:10  kbk
!  possible to store calculated mean fields
!
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0

!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"

   n3d = n3d + 1
   if (n3d .eq. 1) then
      call save_grid_ncdf(ncid,save3d)
   end if

   start(1) = n3d
   edges(1) = 1
   ws(1) = secs
   err = nf_put_vara_real(ncid,time_id,start,edges,ws(1))

   start(1) = 1
   start(2) = 1
   start(3) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = 1

!  Short wave radiation
   call cnv_2d(imin,jmin,imax,jmax,az,swrmean,swr_missing, &
               imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid, swrmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  mean friction velocity
    call cnv_2d(imin,jmin,imax,jmax,az,ustarmean,vel_missing, &
                  imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid,ustarmean_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

!  mean standard deviation of friction velocity
   call cnv_2d(imin,jmin,imax,jmax,az,ustar2mean,vel_missing, &
                  imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid,ustar2mean_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = zlen
   edges(4) = 1

!  layer thickness
   if (hmean_id .gt. 0) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,hmean,h_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf_put_vara_real(ncid,hmean_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10
   end if

!  uumean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,uumean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws)
   err = nf_put_vara_real(ncid, uumean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  vvmean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,vvmean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws)
   err = nf_put_vara_real(ncid, vvmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  wmean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,wmean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws)
   err = nf_put_vara_real(ncid, wmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  salt mean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,Smean,salt_missing, &
               imin,imax,jmin,jmax,0,kmax,ws)
   err = nf_put_vara_real(ncid, saltmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  mean temperature
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,Tmean,temp_missing, &
               imin,imax,jmin,jmax,0,kmax,ws)
   err = nf_put_vara_real(ncid, tempmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

   if (save_mix_analysis) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix3d_S_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf_put_vara_real(ncid, nm3dS_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix3d_T_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf_put_vara_real(ncid, nm3dT_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  phymix3d_S_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf_put_vara_real(ncid, pm3dS_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  phymix3d_T_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf_put_vara_real(ncid, pm3dT_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   

      start(1) = 1
      start(2) = 1
      start(3) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = 1

      call cnv_2d(imin,jmin,imax,jmax,az,nummix2d_S_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, nm2dS_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_2d(imin,jmin,imax,jmax,az,nummix2d_T_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, nm2dT_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_2d(imin,jmin,imax,jmax,az,phymix2d_S_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, pm2dS_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      call cnv_2d(imin,jmin,imax,jmax,az,phymix2d_T_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, pm2dT_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   end if

   err = nf_sync(ncid)
   if (err .NE. NF_NOERR) go to 10

   return

10 FATAL 'save_mean_ncdf: ',nf_strerror(err)
   stop

   return
   end subroutine save_mean_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Adolf Stips and Karsten Bolding (BBH)           !
!-----------------------------------------------------------------------
