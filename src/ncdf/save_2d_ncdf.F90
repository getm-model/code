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
   use variables_2d, only: z,D,U,DU,V,DV,res_u,res_v
   use variables_les, only: AmC_2d
#if USE_BREAKS
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

! average zonal velocity
      if (destag) then
         call to_2d_u(imin,jmin,imax,jmax,az,u,DU,vel_missing,         &
                      imin,jmin,imax,jmax,ws)
      else
         call to_2d_vel(imin,jmin,imax,jmax,au,u,DU,vel_missing,       &
                        imin,jmin,imax,jmax,ws)
      endif
      err = nf90_put_var(ncid,u_id,ws(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

! average meridional velocity
      if (destag) then
         call to_2d_v(imin,jmin,imax,jmax,az,v,DV,vel_missing,         &
                      imin,jmin,imax,jmax,ws)
      else
         call to_2d_vel(imin,jmin,imax,jmax,av,v,DV,vel_missing,       &
              imin,jmin,imax,jmax,ws)
      endif
      err = nf90_put_var(ncid,v_id,ws(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

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
      if (Am_method.eq.AM_LES .and. save_Am_2d) then
         call cnv_2d(imin,jmin,imax,jmax,az,AmC_2d,Am_2d_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf90_put_var(ncid,Am_2d_id,ws(_2D_W_),start,edges)
         if (err .NE. NF90_NOERR) go to 10
      end if

   else ! residual velocities

      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen

      call cnv_2d(imin,jmin,imax,jmax,az,res_u,vel_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf90_put_var(ncid,res_u_id,ws(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

      call cnv_2d(imin,jmin,imax,jmax,az,res_v,vel_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf90_put_var(ncid,res_v_id,ws(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10

#if USE_BREAKS
      err = nf90_put_var(ncid,break_stat_id, &
                         break_stat(_2D_W_),start,edges)
      if (err .NE. NF90_NOERR) go to 10
#endif
   end if
   err = nf90_sync(ncid)
   if (err .NE. NF90_NOERR) go to 10

   return

10 FATAL 'save_2d_ncdf: ',nf90_strerror(err)
   stop

   return
   end subroutine save_2d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
