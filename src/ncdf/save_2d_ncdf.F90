!$Id: save_2d_ncdf.F90,v 1.6 2005-04-25 09:32:34 kbk Exp $
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
   use exceptions
   use ncdf_2d
   use grid_ncdf
   use domain,       only: ioff,joff,imin,imax,jmin,jmax
   use domain,       only: H,az,au,av,min_depth
   use variables_2d, only: z,D,U,DU,V,DV,res_u,res_v,surfdiv
   use meteo,        only: metforcing,calc_met
   use meteo,        only: airp,u10,v10,t2,hum,tcc
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
!  $Log: save_2d_ncdf.F90,v $
!  Revision 1.6  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.5  2003/06/17 14:53:29  kbk
!  default meteo variables names comply with Adolf Stips suggestion + southpole(3)
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
!  Revision 1.3  2001/10/26 12:18:06  bbh
!  No actual storing of data in init_2d_ncdf.F90 -> save_2d_ncdf.F90
!
!  Revision 1.2  2001/09/27 08:35:10  bbh
!  Saving meteo again - in .2d.nc file
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: start(3),edges(3)
   integer, save             :: n2d=0

!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"

   if (secs .ge. _ZERO_) then

      n2d = n2d + 1
      if (n2d .eq. 1) then
         call save_grid_ncdf(ncid,save3d,x_dim,y_dim)
      end if

      start(1) = n2d
      edges(1) = 1
      ws(1) = secs
      err = nf_put_vara_real(ncid,time_id,start,edges,ws(1))

      start(1) = 1
      start(2) = 1
      start(3) = n2d
      edges(1) = xlen
      edges(2) = ylen
      edges(3) = 1

! elevations
      call eta_mask(imin,jmin,imax,jmax,az,H,D,z, &
                    min_depth,elev_missing,imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid,elev_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

! average zonal velocity
      if (destag) then
         call to_2d_u(imin,jmin,imax,jmax,az,u,DU,vel_missing,         &
                      imin,jmin,imax,jmax,ws)
      else
         call to_2d_vel(imin,jmin,imax,jmax,au,u,DU,vel_missing,       &
                        imin,jmin,imax,jmax,ws)
      endif
      err = nf_put_vara_real(ncid, u_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

! average meridional velocity
      if (destag) then
         call to_2d_v(imin,jmin,imax,jmax,az,v,DV,vel_missing,         &
                      imin,jmin,imax,jmax,ws)
      else
         call to_2d_vel(imin,jmin,imax,jmax,av,v,DV,vel_missing,       &
              imin,jmin,imax,jmax,ws)
      endif
      err = nf_put_vara_real(ncid, v_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

! divergence
      call cnv_2d(imin,jmin,imax,jmax,az,surfdiv,divergence_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, surfdiv_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

      if (metforcing .and. save_meteo) then

         if (calc_met) then
            call cnv_2d(imin,jmin,imax,jmax,az,u10,vel_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, u10_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,v10,vel_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, v10_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,airp,airp_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, airp_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,t2,t2_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, t2_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,hum,hum_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, hum_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(imin,jmin,imax,jmax,az,tcc,tcc_missing, &
                        imin,jmin,imax,jmax,ws)
            err = nf_put_vara_real(ncid, tcc_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10
         end if

         call cnv_2d(imin,jmin,imax,jmax,az,tausx,stress_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf_put_vara_real(ncid, tausx_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,tausy,stress_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf_put_vara_real(ncid, tausy_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,swr,swr_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf_put_vara_real(ncid, swr_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(imin,jmin,imax,jmax,az,shf,shf_missing, &
                     imin,jmin,imax,jmax,ws)
         err = nf_put_vara_real(ncid, shf_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

      end if

   else ! residual velocities

      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen

      call cnv_2d(imin,jmin,imax,jmax,az,res_u,vel_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, res_u_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

      call cnv_2d(imin,jmin,imax,jmax,az,res_v,vel_missing, &
                  imin,jmin,imax,jmax,ws)
      err = nf_put_vara_real(ncid, res_v_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   end if
   err = nf_sync(ncid)
   if (err .NE. NF_NOERR) go to 10

   return

10 FATAL 'save_2d_ncdf: ',nf_strerror(err)
   stop

   return
   end subroutine save_2d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
