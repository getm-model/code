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
   use ncdf_3d
   use grid_ncdf,    only: xlen,ylen,zlen
   use domain,       only: ioff,joff,imin,imax,jmin,jmax,kmax
   use domain,       only: H,HU,HV,az,au,av,min_depth
   use domain,       only: convc
#if defined CURVILINEAR || defined SPHERICAL
   use domain,       only: dxv,dyu,arcd1
#else
   use domain,       only: dx,dy,ard1
#endif
   use variables_2d, only: z,D
   use variables_2d, only: U,V,DU,DV
   use variables_3d, only: dt,kmin,ho,hn,uu,hun,vv,hvn,ww,hcc,SS
   use variables_3d, only: taubx,tauby
#ifdef _MOMENTUM_TERMS_
   use variables_3d, only: tdv_u,adv_u,vsd_u,hsd_u,cor_u,epg_u,ipg_u
   use variables_3d, only: tdv_v,adv_v,vsd_v,hsd_v,cor_v,epg_v,ipg_v
#endif
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho,rad,NN
#endif
   use variables_3d, only: nummix3d_S,nummix3d_T,phymix3d_S,phymix3d_T
   use variables_3d, only: numdis3d
   use variables_3d, only: tke,num,nuh,eps
#ifdef SPM
   use variables_3d, only: spm_pool,spm
#endif
#ifdef SPM
   use suspended_matter, only: spm_save
#endif
#ifdef GETM_BIO
   use bio_var, only: numc
   use variables_3d, only: cc3d
#endif
#ifdef _FABM_
   use getm_fabm,only: model,fabm_pel,fabm_ben,fabm_diag,fabm_diag_hz
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
! !LOCAL VARIABLES:
   integer                   :: err,n
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0
   REALTYPE                  :: DONE(E2DFIELD)
   REALTYPE                  :: dum(1)
   integer                   :: i,j
   REALTYPE                  :: uutmp(I3DFIELD),vvtmp(I3DFIELD)
#if defined(CURVILINEAR)
   REALTYPE                  :: uurot(I3DFIELD),vvrot(I3DFIELD)
   REALTYPE                  :: deg2rad = 3.141592654/180.
   REALTYPE                  :: cosconv,sinconv
#endif
   REALTYPE,dimension(E2DFIELD) :: ws2d
   REALTYPE,dimension(I3DFIELD) :: ws
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
      err = nf90_put_var(ncid,hcc_id,ws(_3D_W_),start(1:3),edges(1:3))
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
                       hun,uu,vel_missing,uutmp)
         call to_3d_vv (imin,jmin,imax,jmax,kmin,kmax,az, &
                        hvn,vv,vel_missing,vvtmp)
      else
         call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,au, &
                        hun,uu,vel_missing,uutmp)
         call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,av, &
                        hvn,vv,vel_missing,vvtmp)
      endif
      err = nf90_put_var(ncid,uu_id,uutmp(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_put_var(ncid,vv_id,vvtmp(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      call tow(imin,jmin,imax,jmax,kmin,kmax,az, &
               dt,                               &
#if defined CURVILINEAR || defined SPHERICAL
               dxv,dyu,arcd1,                    &
#else
               dx,dy,ard1,                       &
#endif
               H,HU,HV,hn,ho,uu,hun,vv,hvn,ww,vel_missing,ws)
      err = nf90_put_var(ncid,w_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

#ifdef _MOMENTUM_TERMS_
      err = nf90_put_var(ncid,tdv_u_id,tdv_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,adv_u_id,adv_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,vsd_u_id,vsd_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,hsd_u_id,hsd_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,cor_u_id,cor_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,epg_u_id,epg_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,ipg_u_id,ipg_u(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,tdv_v_id,tdv_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,adv_v_id,adv_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,vsd_v_id,vsd_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,hsd_v_id,hsd_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,cor_v_id,cor_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,epg_v_id,epg_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_put_var(ncid,ipg_v_id,ipg_v(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
#endif

#if defined(CURVILINEAR)
! rotated zonal and meridional velocities
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .gt. 0) then
               cosconv = cos(-convc(i,j)*deg2rad)
               sinconv = sin(-convc(i,j)*deg2rad)
               uurot(i,j,:) = uutmp(i,j,:)*cosconv-vvtmp(i,j,:)*sinconv
               vvrot(i,j,:) = uutmp(i,j,:)*sinconv+vvtmp(i,j,:)*cosconv
            else
               uurot(i,j,:) = vel_missing
               vvrot(i,j,:) = vel_missing
            end if
         end do
      end do
      err = nf90_put_var(ncid,uurot_id,uurot(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_put_var(ncid,vvrot_id,vvrot(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
#endif

   end if

#ifndef NO_BAROCLINIC

   if (salt_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,S,salt_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,salt_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (temp_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,T,temp_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,temp_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (sigma_t_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,rho-1000.,rho_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,sigma_t_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (rad_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,rad,rad_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,rad_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

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

   if (save_numerical_analyses) then

      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,numdis3d,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws)
      err = nf90_put_var(ncid,nm3d_id,ws(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

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

   end if ! save_numerical_analyses

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

#ifdef _FABM_
    if (allocated(fabm_pel)) then
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = zlen
      edges(4) = 1
      do n=1,size(model%state_variables)
         if (fabm_ids(n)==-1) cycle
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,fabm_pel(:,:,:,n), &
                     model%state_variables(n)%missing_value,imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,fabm_ids(n),ws(_3D_W_),start,edges)
         if (err .NE.  NF90_NOERR) go to 10
      end do
      do n=1,size(model%diagnostic_variables)
         if (fabm_ids_diag(n)==-1) cycle
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,fabm_diag(:,:,:,n), &
                     model%diagnostic_variables(n)%missing_value,imin,imax,jmin,jmax,0,kmax,ws)
         err = nf90_put_var(ncid,fabm_ids_diag(n),ws(_3D_W_),start,edges)
         if (err .NE.  NF90_NOERR) go to 10
      end do
      start(3) = n3d
      edges(3) = 1
      do n=1,size(model%bottom_state_variables)
         if (fabm_ids_ben(n)==-1) cycle
         call cnv_2d(imin,jmin,imax,jmax,az,fabm_ben(:,:,n), &
                     model%bottom_state_variables(n)%missing_value,imin,jmin,imax,jmax,ws2d)
         err = nf90_put_var(ncid,fabm_ids_ben(n),ws2d(_2D_W_),start(1:3),edges(1:3))
         if (err .NE.  NF90_NOERR) go to 10
      end do
      do n=1,size(model%horizontal_diagnostic_variables)
         if (fabm_ids_diag_hz(n)==-1) cycle
         call cnv_2d(imin,jmin,imax,jmax,az,fabm_diag_hz(:,:,n), &
                     model%horizontal_diagnostic_variables(n)%missing_value,imin,jmin,imax,jmax,ws2d)
         err = nf90_put_var(ncid,fabm_ids_diag_hz(n),ws2d(_2D_W_),start(1:3),edges(1:3))
         if (err .NE.  NF90_NOERR) go to 10
      end do
   end if
#endif

   if (sync_3d .ne. 0 .and. mod(n3d,sync_3d) .eq. 0) then
      err = nf90_sync(ncid)
      if (err .NE. NF90_NOERR) go to 10
   end if

   return

10 FATAL 'save_3d_ncdf: ',nf90_strerror(err)
   stop

   return
   end subroutine save_3d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
