!$Id: save_3d_ncdf.F90,v 1.4 2003-05-09 11:38:26 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: save_3d_ncdf() - saves 2D-fields.
!
! !INTERFACE:
   subroutine save_3d_ncdf(secs)
!
! !DESCRIPTION:
!
! !USES:
   use ncdf_3d
   use domain, only: ioff,joff,imin,imax,jmin,jmax
   use domain, only: H,az,au,av,min_depth
#if defined(SPHERICAL)
   use domain, only: lonc,latc
#endif
#if ! defined(SPHERICAL)
   use domain, only: xc,yc
#endif
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use domain, only: grid_type,vert_cord,ga
   use variables_2d, only: z,D,u,DU,v,DV
   use variables_3d, only: kmin,hn,uu,hun,vv,hvn,ww
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho
#endif
   use variables_3d, only: tke,num,nuh,eps
#ifndef NO_SUSP_MATTER
   use variables_3d, only: spm
#endif
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
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: save_3d_ncdf.F90,v $
!  Revision 1.4  2003-05-09 11:38:26  kbk
!  added proper undef support - based on Adolf Stips patch
!
!  Revision 1.3  2003/04/23 11:53:24  kbk
!  save lat/lon info for spherical grid
!
!  Revision 1.2  2003/04/07 12:43:12  kbk
!  SPHERICAL and NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:01:48  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/10/25 16:16:21  bbh
!  No actual storing of data in init_3d_ncdf.F90 -> save_3d_ncdf.F90
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
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0
   REALTYPE, parameter       :: x=-1025./9.82
   integer                   :: i,j,k,itmp(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"
   n3d = n3d + 1
   if (n3d .eq. 1) then

      if( xlen*ylen*zlen .gt. size_3d ) then
         FATAL 'Increase size_3d in ncdf_3d() - this needs a fix'
         stop 'ncdf_3d_save'
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

!  elevations
   call eta_mask(imin,jmin,imax,jmax,az,H,D,z,min_depth,elev_missing, &
                 iimin,jmin,iimax,jjmax,ws)
   err = nf_put_vara_real(ncid,elev_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

!  depth integrated zonal velocity
   call to_2d_vel(imin,jmin,imax,jmax,au,u,DU,vel_missing, &
                  imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid,u_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

!  depth integrated meridional velocity
   call to_2d_vel(imin,jmin,imax,jmax,av,v,DV,vel_missing, &
                  imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid,v_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = zlen
   edges(4) = 1

   if (h_id .gt. 0) then
      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                  kmin,az,hn,h_missing,ws)
      err = nf_put_vara_real(ncid,h_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10
   end if

   if (save_vel) then
      call to_3d_vel(imin,jmin,imax,jmax,au, &
                     iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,hn,uu,vel_missing,ws)
      err = nf_put_vara_real(ncid,uu_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

      call to_3d_vel(imin,jmin,imax,jmax,av, &
                     iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,hn,vv,vel_missing,ws)
      err = nf_put_vara_real(ncid,vv_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

      call tow(ws,ww,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
      err = nf_put_vara_real(ncid,w_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10
   end if

#ifndef NO_BAROCLINIC
   if (save_strho) then

      if (save_s) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,S,salt_missing,ws)
         err = nf_put_vara_real(ncid, salt_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_t) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,T,temp_missing,ws)
         err = nf_put_vara_real(ncid, temp_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_rho) then
!KBK        call cnv_3d(ws,x*rho+drho,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
!KBK         call cnv_3d(ws,x*rho+25.,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,x*rho+25.,rho_missing,ws)
         err = nf_put_vara_real(ncid, sigma_t_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if
   end if ! save_strho
#endif

   if (save_turb) then

      if (save_tke) then
!KBK         call cnv_3d(ws,tke,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,tke,tke_missing,ws)
         err = nf_put_vara_real(ncid,tke_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_num) then
!KBK         call cnv_3d(ws,num,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,num,num_missing,ws)
         err = nf_put_vara_real(ncid,num_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_nuh) then
!KBK         call cnv_3d(ws,nuh,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,nuh,nuh_missing,ws)
         err = nf_put_vara_real(ncid,nuh_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_eps) then
!KBK         call cnv_3d(ws,eps,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,eps,eps_missing,ws)
         err = nf_put_vara_real(ncid, eps_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if
   end if ! save_turb

#ifndef NO_BAROCLINIC
   if (save_spm) then
!KBK      call cnv_3d(ws,spm,iimin,jjmin,0,iimax,jjmax,kmax,size_3d)
      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                  kmin,az,spm,spm_missing,ws)
      err = nf_put_vara_real(ncid, spm_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   end if
#endif

   err = nf_sync(ncid)
   if (err .NE. NF_NOERR) go to 10

   return

10 FATAL 'ncdf_3d_save: ',nf_strerror(err)
   stop

   return
   end subroutine save_3d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
