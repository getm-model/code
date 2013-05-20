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
   use netcdf
   use exceptions
   use grid_ncdf,    only: xlen,ylen,zlen,h_missing
   use ncdf_2d, only: ws2d => ws
   use ncdf_3d, only: ws3d => ws
   use ncdf_mean
   use diagnostic_variables
   use domain,       only: ioff,joff,imin,imax,jmin,jmax,kmax
   use domain,       only: H,az,au,av
   use variables_3d, only: kmin
   use m3d, only: calc_temp,calc_salt
#ifdef GETM_BIO
   use bio_var, only: numc
#endif
#ifdef _FABM_
   use gotm_fabm, only: model
#endif

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
! !LOCAL VARIABLES:
   integer                   :: n
   integer                   :: err
   integer                   :: start(4),edges(4)
   integer, save             :: n3d=0
   REALTYPE                  :: dum(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   n3d = n3d + 1
   if (n3d .eq. 1) then
      call save_grid_ncdf(ncid,save3d)
   end if

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

!  Short wave radiation
   call cnv_2d(imin,jmin,imax,jmax,az,swrmean,swr_missing, &
               imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,swrmean_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  mean friction velocity
    call cnv_2d(imin,jmin,imax,jmax,az,ustarmean,vel_missing, &
                  imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,ustarmean_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  mean standard deviation of friction velocity
   call cnv_2d(imin,jmin,imax,jmax,az,ustar2mean,vel_missing, &
                  imin,jmin,imax,jmax,ws2d)
   err = nf90_put_var(ncid,ustar2mean_id,ws2d(_2D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

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
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,hmean_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

!  uumean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,uumean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws3d)
   err = nf90_put_var(ncid, uumean_id,ws3d(_3D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  vvmean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,vvmean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws3d)
   err = nf90_put_var(ncid, vvmean_id,ws3d(_3D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

!  wmean
   call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,wmean,vel_missing, &
               imin,imax,jmin,jmax,0,kmax,ws3d)
   err = nf90_put_var(ncid, wmean_id,ws3d(_3D_W_),start,edges)
   if (err .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
!  salt mean
   if (saltmean_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,Smean,salt_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid, saltmean_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

!  mean temperature
   if (tempmean_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,Tmean,temp_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid, tempmean_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

!  mean sigma_t
   if (sigma_tmean_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,rhomean-1000.,rho_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,sigma_tmean_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
#endif

   if (nd3d_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  numdis_3d_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nd3d_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nd3do_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  numdis_3d_old_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nd3do_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pd3d_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  phydis_3d_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,pd3d_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmS_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix_S_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nmS_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmSo_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix_S_old_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nmSo_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pmS_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  phymix_S_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,pmS_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmT_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix_T_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nmT_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmTo_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  nummix_T_old_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,nmTo_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pmT_id .ne. -1) then
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az, &
                  phymix_T_mean,nummix_missing, &
                  imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid,pmT_id,ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

   start(1) = 1
   start(2) = 1
   start(3) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = 1

   if (ndint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,numdis_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,ndint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pdint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,phydis_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,pdint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmSint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,nummix_S_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,nmSint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pmSint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,phymix_S_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,pmSint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (nmTint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,nummix_T_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,nmTint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (pmTint_id .ne. -1) then
      call cnv_2d(imin,jmin,imax,jmax,az,phymix_T_int_mean,nummix_missing, &
                  imin,jmin,imax,jmax,ws2d)
      err = nf90_put_var(ncid,pmTint_id,ws2d(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end if

#ifdef GETM_BIO
   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = n3d
   edges(1) = xlen
   edges(2) = ylen
   edges(3) = zlen
   edges(4) = 1

   do n=1,numc
      call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,cc3dmean(n,:,:,:), &
                  bio_missing,imin,imax,jmin,jmax,0,kmax,ws3d)
      err = nf90_put_var(ncid, biomean_id(n), ws3d(_3D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
   end do
#endif
#ifdef _FABM_
    if (allocated(fabmmean_pel)) then
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = n3d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = zlen
      edges(4) = 1
      do n=1,size(model%info%state_variables)
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,fabmmean_pel(:,:,:,n), &
                     model%info%state_variables(n)%missing_value,imin,imax,jmin,jmax,0,kmax,ws3d)
         err = nf90_put_var(ncid,fabmmean_ids(n),ws3d(_3D_W_),start,edges)
         if (err .NE.  NF90_NOERR) go to 10
      end do
      do n=1,size(model%info%diagnostic_variables)
         call cnv_3d(imin,jmin,imax,jmax,kmin,kmax,az,fabmmean_diag(:,:,:,n), &
                     model%info%diagnostic_variables(n)%missing_value,imin,imax,jmin,jmax,0,kmax,ws3d)
         err = nf90_put_var(ncid,fabmmean_ids_diag(n),ws3d(_3D_W_),start,edges)
         if (err .NE.  NF90_NOERR) go to 10
      end do
      start(3) = n3d
      edges(3) = 1
      do n=1,size(model%info%state_variables_ben)
         call cnv_2d(imin,jmin,imax,jmax,az,fabmmean_ben(:,:,n), &
                     model%info%state_variables_ben(n)%missing_value,imin,jmin,imax,jmax,ws2d)
         err = nf90_put_var(ncid,fabmmean_ids_ben(n),ws2d(_2D_W_),start(1:3),edges(1:3))
         if (err .NE.  NF90_NOERR) go to 10
      end do
      do n=1,size(model%info%diagnostic_variables_hz)
         call cnv_2d(imin,jmin,imax,jmax,az,fabmmean_diag_hz(:,:,n), &
                     model%info%diagnostic_variables_hz(n)%missing_value,imin,jmin,imax,jmax,ws2d)
         err = nf90_put_var(ncid,fabmmean_ids_diag_hz(n),ws2d(_2D_W_),start(1:3),edges(1:3))
         if (err .NE.  NF90_NOERR) go to 10
      end do
   end if
#endif

   err = nf90_sync(ncid)
   if (err .NE. NF90_NOERR) go to 10

   return

10 FATAL 'save_mean_ncdf: ',nf90_strerror(err)
   stop

   return
   end subroutine save_mean_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Adolf Stips and Karsten Bolding (BBH)           !
!-----------------------------------------------------------------------
