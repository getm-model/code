#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: save_2d_ncdf() - saves 2D-fields.
!
! !INTERFACE:
   subroutine save_2d_ncdf(secs)
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use exceptions
   use ncdf_2d
   use grid_ncdf,    only: xlen,ylen
   use domain,       only: ioff,joff,imin,imax,jmin,jmax
   use domain,       only: H,az,au,av,crit_depth
   use domain,       only: convc
   use variables_2d, only: z,D,U,DU,V,DV,res_u,res_v
#ifdef USE_BREAKS
   use variables_2d, only: break_stat
#endif
   use meteo,        only: metforcing,calc_met
   use meteo,        only: airp,u10,v10,t2,hum,tcc
   use meteo,        only: evap,precip
   use meteo,        only: tausx,tausy,swr,shf

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: secs
!
! !DEFINED PARAMTERS:
   logical, parameter                  :: save3d=.false.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: start(3),edges(3)
   integer, save             :: n2d=0
   REALTYPE                  :: dum(1)
   integer                   :: i,j
   REALTYPE                  :: Utmp(E2DFIELD),Vtmp(E2DFIELD)
#if defined(CURVILINEAR)
   REALTYPE                  :: Urot(E2DFIELD),Vrot(E2DFIELD)
   REALTYPE                  :: deg2rad = 3.141592654/180.
   REALTYPE                  :: cosconv,sinconv
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
   if (secs .ge. _ZERO_) then

      n2d = n2d + 1
      if (n2d .eq. 1) then
         call save_grid_ncdf(ncid,save3d)
      end if

      start(1) = n2d
      edges(1) = 1
      dum(1) = secs
      err = nf90_put_var(ncid,time_id,dum,start,edges)
      if (err .NE. NF90_NOERR) go to 10

      start(1) = 1
      start(2) = 1
      start(3) = n2d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = 1

! elevations
      call eta_mask(imin,jmin,imax,jmax,az,H,D,z, &
                    crit_depth,elev_missing,imin,jmin,imax,jmax,ws)
      err = nf90_put_var(ncid,elev_id,ws(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

! grid-related zonal and meridional velocities
      if (destag) then
         call to_2d_u(imin,jmin,imax,jmax,az,U,DU,vel_missing,      &
                      imin,jmin,imax,jmax,Utmp)
         call to_2d_v(imin,jmin,imax,jmax,az,V,DV,vel_missing,      &
                      imin,jmin,imax,jmax,Vtmp)
      else
         call to_2d_vel(imin,jmin,imax,jmax,au,U,DU,vel_missing,       &
                        imin,jmin,imax,jmax,Utmp)
         call to_2d_vel(imin,jmin,imax,jmax,av,V,DV,vel_missing,       &
                        imin,jmin,imax,jmax,Vtmp)
      endif
      err = nf90_put_var(ncid,u_id,Utmp(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_put_var(ncid,v_id,Vtmp(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

#if defined(CURVILINEAR)
! rotated zonal and meridional velocities
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .gt. 0) then
               cosconv = cos(-convc(i,j)*deg2rad)
               sinconv = sin(-convc(i,j)*deg2rad)
               Urot(i,j) = Utmp(i,j)*cosconv-Vtmp(i,j)*sinconv
               Vrot(i,j) = Utmp(i,j)*sinconv+Vtmp(i,j)*cosconv
            else
               Urot(i,j) = vel_missing
               Vrot(i,j) = vel_missing
            end if
         end do
      end do
      err = nf90_put_var(ncid,urot_id,Urot(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_put_var(ncid,vrot_id,Vrot(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
#endif

      if (metforcing .and. save_meteo) then

         if (calc_met) then
            call cnv_2d(imin,jmin,imax,jmax,az,u10,vel_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,u10_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,v10,vel_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,v10_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,airp,airp_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,airp_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,t2,t2_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,t2_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,hum,hum_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,hum_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,tcc,tcc_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf90_put_var(ncid,tcc_id,ws(_2D_W_),start,edges)
            if (err .NE. NF90_NOERR) go to 10

            if (evap_id .ge. 0) then
               call cnv_2d(imin,jmin,imax,jmax,az,evap,evap_missing, &
                          imin,jmin,imax,jmax,ws)
               err = nf90_put_var(ncid,evap_id,ws(_2D_W_),start,edges)
               if (err .NE. NF90_NOERR) go to 10
            end if

            if (precip_id .ge. 0) then
               call cnv_2d(imin,jmin,imax,jmax,az,precip,precip_missing, &
                          imin,jmin,imax,jmax,ws)
               err = nf90_put_var(ncid,precip_id,ws(_2D_W_),start,edges)
               if (err .NE. NF90_NOERR) go to 10
            end if

         end if

         call cnv_2d(imin,jmin,imax,jmax,az,tausx,stress_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,tausx_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,tausy,stress_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,tausy_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,swr,swr_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,swr_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,shf,shf_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,shf_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10

      end if

   else ! residual velocities

!     Note (KK): there are conceptual discrepancies in the implementation
!                of the residual transports. therefore the buggy output
!                into 2d ncdf is not fixed (either add missing start(3)
!                and edges(3), or define the ncdf fields independent on
!                time) and the activation with residual.gt.0 is not
!                recommended.

      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen

      if (res_u_id .ne. -1) then
         call cnv_2d(imin,jmin,imax,jmax,az,res_u,vel_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,res_u_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

      if (res_v_id .ne. -1) then
         call cnv_2d(imin,jmin,imax,jmax,az,res_v,vel_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,res_v_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

#ifdef USE_BREAKS
      err = nf90_put_var(ncid,break_stat_id,break_stat(_2D_W_))
      if (err .NE. NF90_NOERR) go to 10
#endif
   end if
   if (sync_2d .ne. 0 .and. mod(n2d,sync_2d) .eq. 0) then
      err = nf90_sync(ncid)
      if (err .NE. NF90_NOERR) go to 10
   end if

   return

10 FATAL 'save_2d_ncdf: ',nf90_strerror(err)
   stop

   return
   end subroutine save_2d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
