!$Id: save_2d_ncdf.F90,v 1.1 2002-05-02 14:01:48 gotm Exp $
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
   use ncdf_2d
   use domain, only: ioff,joff,imin,imax,jmin,jmax,H,au,av,min_depth
   use domain, only: grid_type,xc,yc
   use variables_2d, only: z,D,U,DU,V,DV,res_u,res_v,surfdiv
   use meteo, only: calc_met
   use meteo, only: airp,u10,v10,t2,hum,cc
   use meteo, only: tausx,tausy,swr,shf
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
!  $Log: save_2d_ncdf.F90,v $
!  Revision 1.1  2002-05-02 14:01:48  gotm
!  Initial revision
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
   integer      :: err
   integer	:: start(3),edges(3)
   integer, save        :: n2d=0
   integer	:: i,j,itmp(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"
   if (secs .ge. _ZERO_) then
      n2d = n2d + 1
      if (n2d .eq. 1) then

         if( xlen*ylen .gt. size_2d ) then
            FATAL 'Increase size_2d in ncdf_2d_save() - this needs a fix'
	    stop 'ncdf_2d_save'
         end if

!     save info on offset, grid type and vertical coordinates
         itmp(1) = grid_type
         err = nf_put_var_int(ncid,grid_type_id,itmp)
         if (err .NE. NF_NOERR) go to 10

         itmp(1) = ioff
         err = nf_put_var_int(ncid,ioff_id,itmp)
         if (err .NE. NF_NOERR) go to 10

         itmp(1) = joff
         err = nf_put_var_int(ncid,joff_id,itmp)
         if (err .NE. NF_NOERR) go to 10

!  save coordinate information
         select case (grid_type)
#ifndef CURVILINEAR
            case (1,2)
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
#else
            case (3)
               STDERR 'xc and yc are read from input file directly'
#endif
            case default
         end select

         start(1) = 1
         start(2) = 1
         edges(1) = xlen
         edges(2) = ylen

         call cnv_2d(ws,imin,jmin,imax,jmax,H,imin,jmin,imax,jmax)
         err = nf_put_vara_real(ncid, bathymetry_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         err = nf_sync(ncid)
         if (err .NE. NF_NOERR) go to 10

      end if ! (n2d .eq. 1)

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
      call eta_mask(ws,imin,jmin,imax,jmax,z,D,H,imin,jmin,imax,jmax,min_depth)
      err = nf_put_vara_real(ncid,elev_id,start,edges,ws)
      if (err .NE. NF_NOERR) go to 10

! average zonal velocity
      call to_2d_vel(ws,imin,jmin,imax,jmax,au,u,DU,imin,jmin,imax,jmax)
      err = nf_put_vara_real(ncid, u_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

! average meridional velocity
      call to_2d_vel(ws,imin,jmin,imax,jmax,av,v,DV,imin,jmin,imax,jmax)
      err = nf_put_vara_real(ncid, v_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

! divergence
      call cnv_2d(ws,imin,jmin,imax,jmax,surfdiv,imin,jmin,imax,jmax)
      err = nf_put_vara_real(ncid, surfdiv_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

      if (save_meteo) then

         if (calc_met) then
            call cnv_2d(ws,imin,jmin,imax,jmax,u10,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, u10_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(ws,imin,jmin,imax,jmax,v10,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, v10_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(ws,imin,jmin,imax,jmax,airp,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, airp_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(ws,imin,jmin,imax,jmax,t2,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, t2_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(ws,imin,jmin,imax,jmax,hum,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, hum_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10

            call cnv_2d(ws,imin,jmin,imax,jmax,cc,imin,jmin,imax,jmax)
            err = nf_put_vara_real(ncid, cc_id, start, edges, ws)
            if (err .NE. NF_NOERR) go to 10
         end if

         call cnv_2d(ws,imin,jmin,imax,jmax,tausx,imin,jmin,imax,jmax)
         err = nf_put_vara_real(ncid, tausx_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(ws,imin,jmin,imax,jmax,tausy,imin,jmin,imax,jmax)
         err = nf_put_vara_real(ncid, tausy_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(ws,imin,jmin,imax,jmax,swr,imin,jmin,imax,jmax)
         err = nf_put_vara_real(ncid, swr_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

         call cnv_2d(ws,imin,jmin,imax,jmax,shf,imin,jmin,imax,jmax)
         err = nf_put_vara_real(ncid, shf_id, start, edges, ws)
         if (err .NE. NF_NOERR) go to 10

      end if
   else ! residual velocities
      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen

      call cnv_2d(ws,imin,jmin,imax,jmax,res_u,imin,jmin,imax,jmax)
      err = nf_put_vara_real(ncid, res_u_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10

      call cnv_2d(ws,imin,jmin,imax,jmax,res_v,imin,jmin,imax,jmax)
      err = nf_put_vara_real(ncid, res_v_id, start, edges, ws)
      if (err .NE. NF_NOERR) go to 10
   end if
   err = nf_sync(ncid)
   if (err .NE. NF_NOERR) go to 10

   return

10 FATAL 'ncdf_2d_save: ',nf_strerror(err)
   stop

   return
   end subroutine save_2d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
