!$Id: save_3d_ncdf.F90,v 1.13 2007-02-20 13:52:15 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Save 3D netCDF variables
!
! !INTERFACE:
   subroutine save_3d_ncdf(secs)
!
! !DESCRIPTION:
!
! !USES:
   use exceptions
   use ncdf_3d
   use grid_ncdf
   use domain,       only: ioff,joff,imin,imax,jmin,jmax
   use domain,       only: iimin,iimax,jjmin,jjmax,kmax
   use domain,       only: H,az,au,av,min_depth
   use variables_2d, only: z,D
   use variables_2d, only: U,V,DU,DV
   use variables_3d, only: kmin,hn,uu,hun,vv,hvn,ww,hcc
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho,rad
#endif
   use variables_3d, only: tke,num,nuh,eps
#ifdef SPM
   use variables_3d, only: spm_pool,spm
#endif
#ifdef SPM
   use suspended_matter, only: spm_save
#endif
#ifdef GETM_BIO
   use bio_var, only: numc
   use variables_3d, only: cc3d,ws3d
#endif
   use parameters,   only: g,rho_0
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: secs
!
! !DEFINED PARAMTERS:
   logical, parameter   :: save3d=.true.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: save_3d_ncdf.F90,v $
!  Revision 1.13  2007-02-20 13:52:15  kbk
!  solar radiation -> 3d field - possible to save
!
!  Revision 1.12  2006-03-17 11:06:33  kbk
!  cleaner inclusion of SPM module
!
!  Revision 1.11  2005/09/23 11:27:10  kbk
!  support for biology via GOTMs biology modules
!
!  Revision 1.10  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.9  2004/06/15 08:25:57  kbk
!  added supoort for spm - Ruiz
!
!  Revision 1.8  2004/05/04 09:23:51  kbk
!  hydrostatic consistency criteria stored in .3d.nc file
!
!  Revision 1.7  2003/12/16 12:47:11  kbk
!  rho_0 and g from parameters (manuel)
!
!  Revision 1.6  2003/12/08 07:21:53  hb
!  use proper layer heights for saving velocities
!
!  Revision 1.5  2003/05/09 11:53:13  kbk
!  forgot to delete some debug lines
!
!  Revision 1.4  2003/05/09 11:38:26  kbk
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
   integer                   :: err,n
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0
   REALTYPE, parameter       :: x=-rho_0/g

!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"

   n3d = n3d + 1
   if (n3d .eq. 1) then

      call save_grid_ncdf(ncid,save3d)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = zlen

      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                  kmin,az,hcc,-_ONE_,ws)
      err = nf_put_vara_real(ncid,hcc_id,start,edges,ws)
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
   call to_2d_vel(imin,jmin,imax,jmax,au,U,DU,vel_missing, &
                  imin,jmin,imax,jmax,ws)
   err = nf_put_vara_real(ncid,u_id,start,edges,ws)
   if (err .NE. NF_NOERR) go to 10

!  depth integrated meridional velocity
   call to_2d_vel(imin,jmin,imax,jmax,av,V,DV,vel_missing, &
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
      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax,    &
                  kmin,az,hn,hh_missing,ws)
      err = nf_put_vara_real(ncid,h_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10
   end if

   if (save_vel) then

      if (destag) then
         call to_3d_uu(imin, jmin, imax, jmax,  az,                    &
                       iimin,jjmin,iimax,jjmax,kmax,                   &
                       kmin,hun,uu,vel_missing,ws)
      else
         call to_3d_vel( imin, jmin, imax, jmax,  au,                  &
                        iimin,jjmin,iimax,jjmax,kmax,                  &
                        kmin,hun,uu,vel_missing,ws)
      endif

      err = nf_put_vara_real(ncid,uu_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10


      if (destag) then
         call to_3d_vv ( imin, jmin, imax, jmax,  az,                 &
                        iimin,jjmin,iimax,jjmax,kmax,                 &
                        kmin,hvn,vv,vel_missing,ws)
      else
         call to_3d_vel( imin, jmin, imax, jmax,  av,                 &
                        iimin,jjmin,iimax,jjmax,kmax,                 &
                        kmin,hvn,vv,vel_missing,ws)
      endif

      err = nf_put_vara_real(ncid,vv_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

      call tow(ws,ww,iimin,jjmin,0,iimax,jjmax,kmax)
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
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,x*rho+rho_0-1000.,rho_missing,ws)
         err = nf_put_vara_real(ncid, sigma_t_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_rad) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,rad,rad_missing,ws)
         err = nf_put_vara_real(ncid, rad_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if

   end if ! save_strho
#endif

   if (save_turb) then

      if (save_tke) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,tke,tke_missing,ws)
         err = nf_put_vara_real(ncid,tke_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_num) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,num,num_missing,ws)
         err = nf_put_vara_real(ncid,num_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_nuh) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,nuh,nuh_missing,ws)
         err = nf_put_vara_real(ncid,nuh_id,start,edges,ws)
         if (err .NE. NF_NOERR) go to 10
      end if

      if (save_eps) then
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,eps,eps_missing,ws)
         err = nf_put_vara_real(ncid, eps_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10
      end if
   end if ! save_turb

#ifdef SPM
   if (spm_save) then
      call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax,    &
                  kmin,az,spm,spm_missing,ws)
      err = nf_put_vara_real(ncid, spm_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
      !spm pool is a 2d magnitude
      start(1) = 1
      start(2) = 1
      start(3) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = 1
      call cnv_2d(imin,jmin,imax,jmax,az,spm_pool,spmpool_missing,    &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, spmpool_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   end if
#endif

#ifdef GETM_BIO
!   if (save_bio) then
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = zlen
      edges(4) = 1
      do n=1,numc
         call cnv_3d(imin,jmin,imax,jmax,iimin,jjmin,iimax,jjmax,kmax, &
                     kmin,az,cc3d(n,:,:,:),bio_missing,ws)
         err = nf_put_vara_real(ncid, bio_ids(n), start, edges, ws)
         if (err .NE.  NF_NOERR) go to 10
      end do
!   end if
#endif

   err = nf_sync(ncid)
   if (err .NE. NF_NOERR) go to 10

   return

10 FATAL 'save_3d_ncdf: ',nf_strerror(err)
   stop

   return
   end subroutine save_3d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
