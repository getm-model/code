!$Id: save_3d_ncdf.F90,v 1.19 2009-04-22 10:09:36 lars Exp $
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
   use netcdf
   use exceptions
   use ncdf_2d,      only: ws2d => ws
   use ncdf_3d
   use grid_ncdf,    only: xlen,ylen,zlen
   use domain,       only: ioff,joff,imin,imax,jmin,jmax,kmax
   use domain,       only: H,HU,HV,az,au,av,min_depth
#if defined CURVILINEAR || defined SPHERICAL
   use domain,       only: dxc,dyc
#else
   use domain,       only: dx,dy
#endif
   use variables_2d, only: z,D
   use variables_2d, only: U,V,DU,DV
   use variables_3d, only: dt,kmin,ho,hn,uu,hun,vv,hvn,ww,hcc,SS
   use variables_3d, only: taubx,tauby
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho,rad,NN
   use variables_3d, only: nummix3d_S,nummix3d_T,phymix3d_S,phymix3d_T
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
   use m3d, only: calc_temp,calc_salt
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
!  Revision 1.19  2009-04-22 10:09:36  lars
!  support for bottom stress output
!
!  Revision 1.18  2009-01-05 09:57:06  kb
!  option for storing SS and NN
!
!  Revision 1.17  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.16  2007-05-26 12:19:31  kbk
!  destag of vertical velocities
!
!  Revision 1.15  2007-05-22 09:37:20  kbk
!  saving physical vertical velocities
!
!  Revision 1.14  2007-02-23 12:20:37  kbk
!  introduce buoy() and rename rho() to buoy()i where appropriate
!
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
   REALTYPE                  :: DONE(E2DFIELD)
   REALTYPE                  :: dum(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   n3d = n3d + 1
   if (n3d .eq. 1) then

      call save_grid_ncdf(ncid,save3d)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = zlen
      edges(4) = 1

      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,hcc,-_ONE_, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,hcc_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_sync(ncid)
      if (err .NE. NF90_NOERR) go to 10

   end if ! (n3d .eq. 1)

   start(1) = n3d
   edges(1) = 1
   dum(1) = secs
   err = nf90_put_var(ncid,time_id,dum,start,edges)

   start(1) = 1
   start(2) = 1
   start(3) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = 1

!  elevations
   call eta_mask(imin,jmin,imax,jmax,az,H,D,z,min_depth,elev_missing, &
                 imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,elev_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  depth integrated zonal velocity
   call to_2d_vel(imin,jmin,imax,jmax,au,U,DU,vel_missing, &
                  imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,u_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  depth integrated meridional velocity
   call to_2d_vel(imin,jmin,imax,jmax,av,V,DV,vel_missing, &
                  imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,v_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

   
   if (save_taub) then

      !  bottom stress (x)
      if (destag) then
         DONE = _ONE_
         call to_2d_u(imin,jmin,imax,jmax,au,rho_0*taubx,DONE,tau_missing, &
              imin,jmin,imax,jmax,ws2d)
      else
         call cnv_2d(imin,jmin,imax,jmax,au,rho_0*taubx,tau_missing,       &
              imin,jmin,imax,jmax,ws2d)

      endif

      err = nf90_put_var(ncid,taubx_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10


      !  bottom stress (y)
      if (destag) then
         DONE = _ONE_
         call to_2d_v(imin,jmin,imax,jmax,av,rho_0*tauby,DONE,tau_missing, &
              imin,jmin,imax,jmax,ws2d)
      else
         call cnv_2d(imin,jmin,imax,jmax,av,rho_0*tauby,tau_missing,       &
              imin,jmin,imax,jmax,ws2d)
      endif

      err = nf90_put_var(ncid,tauby_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10


   endif


   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = zlen
   edges(4) = 1

   if (h_id .gt. 0) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,hn,hh_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,h_id,ws(_3D_W_),start,edges)
!      err = nf90_put_var(ncid,h_id,hn(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (save_vel) then

      if (destag) then
         call to_3d_uu(imin,jmin,imax,jmax,kmin,kmax,az, &
                       hun,uu,vel_missing,ws)
      else
         call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,au, &
                        hun,uu,vel_missing,ws)
      endif

      err = nf90_put_var(ncid,uu_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10


      if (destag) then
         call to_3d_vv (imin,jmin,imax,jmax,kmin,kmax,az, &
                        hvn,vv,vel_missing,ws)
      else
         call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,av, &
                        hvn,vv,vel_missing,ws)
      endif

      err = nf90_put_var(ncid,vv_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      call tow(imin,jmin,imax,jmax,kmin,kmax,az, &
               dt,                               &
#if defined CURVILINEAR || defined SPHERICAL
               dxc,dyc,                          &
#else
               dx,dy,                            &
#endif
               HU,HV,hn,ho,uu,hun,vv,hvn,ww,vel_missing,destag,ws)
      err = nf90_put_var(ncid,w_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

   end if

#ifndef NO_BAROCLINIC
   if (save_strho) then

      if (save_s) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,S,salt_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,salt_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_t) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,T,temp_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,temp_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_rho) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,rho-1000.,rho_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,sigma_t_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_rad) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,rad,rad_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,rad_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

   end if ! save_strho
#endif

   if (save_turb) then

      if (save_tke) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,tke,tke_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,tke_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_num) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,num,num_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,num_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_nuh) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,nuh,nuh_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,nuh_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (save_eps) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,eps,eps_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,eps_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if
   end if ! save_turb

   if (save_ss_nn) then

      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,SS,SS_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,SS_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,NN,NN_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,NN_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

#endif

   end if ! save_ss_nn

#ifndef NO_BAROCLINIC
   if (save_mix_analysis) then
      if (calc_salt) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,nummix3d_S,nummix_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,nm3dS_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,phymix3d_S,nummix_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,pm3dS_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (calc_temp) then
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,nummix3d_T,nummix_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,nm3dT_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,phymix3d_T,nummix_missing, &
                     imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,pm3dT_id,ws(_3D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if
   end if ! save_mix_analysis
#endif

#ifdef SPM
   if (spm_save) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,spm,spm_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,spm_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      !spm pool is a 2d magnitude
      start(1) = 1
      start(2) = 1
      start(3) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = 1
      call cnv_2d(imin,jmin,imax,jmax,az,spm_pool,spmpool_missing,    &
                  imin,jmin,imax,jmax,ws)
      err = nf90_put_var(ncid,spmpool_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
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
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,cc3d(n,:,:,:), &
                     bio_missing,imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,bio_ids(n),ws(_3D_W_),start,edges)
         if (err .NE.  NF90_NOERR) go to 10
      end do
!   end if
#endif

   err = nf90_sync(ncid)
   if (err .NE. NF90_NOERR) go to 10

   return

10 FATAL 'save_3d_ncdf: ',nf90_strerror(err)
   stop

   return
   end subroutine save_3d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
