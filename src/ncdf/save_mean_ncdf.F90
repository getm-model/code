!$Id: save_mean_ncdf.F90,v 1.1 2004-03-29 15:38:10 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: save_mean_ncdf() - saves mean-fields.
!
! !INTERFACE:
   subroutine save_mean_ncdf(secs)
!
! !DESCRIPTION:
!
! !USES:
   use ncdf_mean
   use domain, only: ioff,joff,imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use domain, only: H,az
   use domain, only: grid_type,vert_cord,ga
   use variables_3d, only: kmin
#if defined(SPHERICAL)
   use domain, only: lonc,latc
#endif
#if ! defined(SPHERICAL)
   use domain, only: xc,yc
#endif
   use diagnostic_variables
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: secs
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips & Karsten Bolding
!
!  $Log: save_mean_ncdf.F90,v $
!  Revision 1.1  2004-03-29 15:38:10  kbk
!  possible to store calculated mean fields
!
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0
   integer                   :: i,j,k,itmp(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"
   n3d = n3d + 1
   if (n3d .eq. 1) then

      if( xlen*ylen*zlen .gt. size_3d ) then
         FATAL 'Increase size_3d in ncdf_save_mean() - this needs a fix'
         stop 'ncdf_mean_save'
      end if

!     save info on offset, grid type and vertical coordinates
      itmp(1) = grid_type
      err = nf_put_var_int(ncid,grid_type_id,itmp)
      if (err .NE. NF_NOERR) go to 10

      itmp(1) = vert_cord
      err = nf_put_var_int(ncid,vert_cord_id,itmp)
      if (err .NE. NF_NOERR) go to 10

      itmp(1) = ioff
      err = nf_put_var_int(ncid,ioff_id,itmp)
      if (err .NE. NF_NOERR) go to 10

      itmp(1) = joff
      err = nf_put_var_int(ncid,joff_id,itmp)
      if (err .NE. NF_NOERR) go to 10

!     save coordinate information
      select case (grid_type)
         case (1)
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
            do i=imin,imax
               ws(i) = xc(i)
            end do
            err = nf_put_var_real(ncid,xc_id,ws)
            if (err .NE. NF_NOERR) go to 10
            do j=jmin,jmax
               ws(j) = yc(j)
            end do
            err = nf_put_var_real(ncid,yc_id,ws)
            if (err .NE. NF_NOERR) go to 10
#endif
         case (2)
#if defined(SPHERICAL)
            do i=imin,imax
               ws(i) = lonc(i,1)
            end do
            err = nf_put_var_real(ncid,lonc_id,ws)
            if (err .NE. NF_NOERR) go to 10
            do j=jmin,jmax
               ws(j) = latc(1,j)
            end do
            err = nf_put_var_real(ncid,latc_id,ws)
            if (err .NE. NF_NOERR) go to 10
#endif
         case (3)
#if defined(CURVILINEAR)
            STDERR 'xc and yc are read from input file directly'
#endif
         case default
      end select

      select case (vert_cord)
         case (1,2,3)
            do k=0,kmax
               ws(k+1) = ga(k)
            end do
            err = nf_put_var_real(ncid,z_id,ws)
            if (err .NE. NF_NOERR) go to 10
         case default
      end select

      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen
      call cnv_2d(imin,jmin,imax,jmax,az,H,h_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid,bathymetry_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

      err = nf_sync(ncid)
      if (err .NE. NF_NOERR) go to 10

   end if ! (n3d .eq. 1)

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
      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                  kmin,az,hmean,h_missing,ws)
      err = nf_put_vara_real(ncid,hmean_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10
   end if

!  uumean
   call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
               kmin,az,uumean,vel_missing,ws)
   err = nf_put_vara_real(ncid, uumean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  vvmean
   call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
               kmin,az,vvmean,vel_missing,ws)
   err = nf_put_vara_real(ncid, vvmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  wmean
   call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
               kmin,az,wmean,vel_missing,ws)
   err = nf_put_vara_real(ncid, wmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  salt mean
   call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
               kmin,az,Smean,salt_missing,ws)
   err = nf_put_vara_real(ncid, saltmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

!  mean temperature
   call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
               kmin,az,Tmean,temp_missing,ws)
   err = nf_put_vara_real(ncid, tempmean_id, start, edges, ws)
   if (err .NE. NF_NOERR) go to 10

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
