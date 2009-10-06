!$Id: save_grid_ncdf.F90,v 1.10 2009-10-06 13:11:16 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save grid related variables
!
! !INTERFACE:
   subroutine save_grid_ncdf(ncid,save3d)
!
! !DESCRIPTION:
! This routine saves netCDF variables in an already existing netCDF file
! in save mode with netCDF file-id "{\tt ncid}". The variables saved
! correspond to those GETM variables not changing in time, i.e.\
! grid related variables and bathymetry. If the logical flag
! "{\tt save3d}" evaluates false, no information about the vertical grid
! is saved (e.g.\ if results from a horizontally integrated run are stored).
!
! !USES:
   use exceptions
   use netcdf
   use grid_ncdf
   use domain, only: imin,imax,jmin,jmax
   use domain, only: grid_type,vert_cord
   use domain, only: have_lonlat,have_xy
   use domain, only: ioff,joff
   use domain, only: dx,dy
   use domain, only: dlon,dlat
   use domain, only: xcord,ycord
   use domain, only: xc,yc
   use domain, only: xx,yx
   use domain, only: latc,lonc,convc
   use domain, only: latx,lonx,convx
   use domain, only: latu,latv
   use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,dxx,dyx
!KB   use domain, only: rearth
   use domain, only: H,ga
   use domain, only: az,au,av
   use output, only: save_metrics,save_masks

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)        :: ncid
   logical, intent(in)        :: save3d
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: save_grid_ncdf.F90,v $
!  Revision 1.10  2009-10-06 13:11:16  kb
!  only save bathymetry when az > 0
!
!  Revision 1.9  2009-09-25 12:17:25  kb
!  removed undef [IJ]RANGE
!
!  Revision 1.8  2009-09-23 12:40:00  kb
!  only save latu, latv, convc when have_lonlat=.true. when grid_type=1
!
!  Revision 1.7  2009-09-23 10:11:48  kb
!  rewrite of grid-initialisation, optional grid info saved to file, -DSAVE_HALO, updated documentation
!
!  Revision 1.6  2007-10-16 07:14:35  kbk
!  pseudo coordinate variables for curvi-linear grids
!
!  Revision 1.5  2007-03-30 13:11:00  hb
!  Use of adaptive and hybrid vertical coordinates technically enabled
!
!  Revision 1.4  2007-02-22 08:48:13  kbk
!  possible to save masks (az, au, av)
!
!  Revision 1.3  2006-03-10 08:44:02  kbk
!  fixed saving coordinate variables
!
!  Revision 1.2  2005-11-01 15:44:13  kbk
!  fixed saving of lonc instaed of latc
!
!  Revision 1.1  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   integer                   :: status
   integer                   :: start(2),edges(2)
   integer                   :: id
   character(32)             :: zname
   REALTYPE                  :: ws(E2DFIELD)
!EOP
!------------------------------------------------------------------------
!BOC

! set NetCDF slice information
! set SAVE_HALOS via Makefile
#ifdef SAVE_HALOS
   LEVEL3 'include HALOs in NetCDF hostart files'
   start(1) = 1; edges(1) = (imax+HALO)-(imin-HALO)+1
   start(2) = 1; edges(2) = (jmax+HALO)-(jmin-HALO)+1
#else
   start(1) = 1; edges(1) = imax-imin+1
   start(2) = 1; edges(2) = jmax-jmin+1
#endif

!  get NetCDF variable id's and save corresponding variable

!  grid type
   status = nf90_inq_varid(ncid,'grid_type',id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","grid_type_id -")

   status = nf90_put_var(ncid,id,grid_type)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","grid_type -")

   if (save3d) then
      status = nf90_inq_varid(ncid,'vert_cord',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","vert_cord_id -")

      status = nf90_put_var(ncid,id,vert_cord)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                    "save_grid_ncdf()","vert_cord -")
   endif

!  offset for parallel runs
   status = nf90_inq_varid(ncid,'ioff',id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","ioff_id -")
   status = nf90_put_var(ncid,id,ioff)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","ioff -")

   status = nf90_inq_varid(ncid,'joff',id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","joff_id -")
   status = nf90_put_var(ncid,id,joff)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","joff -")

!  grid-specific settings
   select case (grid_type)

      case (1)

!        grid spacing
         status = nf90_inq_varid(ncid,'dx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dx_id -")
         status = nf90_put_var(ncid,id,dx)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                     "save_grid_ncdf()","dx -")

         status = nf90_inq_varid(ncid,'dy',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dy_id -")
         status = nf90_put_var(ncid,id,dy)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                     "save_grid_ncdf()","dy -")

!        coordinate variables
         status = nf90_inq_varid(ncid,'xc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xc_id -")
         status = nf90_put_var(ncid,id,xcord(imin:imax))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xc -")

         status = nf90_inq_varid(ncid,'yc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yc_id -")
         status = nf90_put_var(ncid,id,ycord(jmin:jmax))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","yc -")

         if ( have_lonlat ) then

            status = nf90_inq_varid(ncid,'lonc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","lonc -")
            status = nf90_put_var(ncid,id,lonc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","lonc -")

            status = nf90_inq_varid(ncid,'latc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latc -")
            status = nf90_put_var(ncid,id,latc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latc -")

            status = nf90_inq_varid(ncid,'convc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","convc -")
            status = nf90_put_var(ncid,id,convc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","convc -")

         end if

      case (2)

!        grid spacing
         status = nf90_inq_varid(ncid,'dlon',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dlon_id -")
         status = nf90_put_var(ncid,id,dlon)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","dlon -")

         status = nf90_inq_varid(ncid,'dlat',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dlat_id -")
         status = nf90_put_var(ncid,id,dlat)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","dlat -")

!        coordinate variables
         status = nf90_inq_varid(ncid,'lonc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lonc_id -")
         status = nf90_put_var(ncid,id,xcord(IRANGE))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","lonc -")

         status = nf90_inq_varid(ncid,'latc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","latc_id -")
         status = nf90_put_var(ncid,id,ycord(JRANGE))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","latc -")
         if ( have_xy ) then

            status = nf90_inq_varid(ncid,'xc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","xc -")
            status = nf90_put_var(ncid,id,xc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","xc -")

            status = nf90_inq_varid(ncid,'yc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","yc -")
            status = nf90_put_var(ncid,id,yc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","yc -")

         end if

#if 1
      case (3,4)

!        pseudo coordinate variables
         status = nf90_inq_varid(ncid,'xic',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xic_id -")
         do i=imin,imax
            xcord(i) = ioff+i
         end do
         status = nf90_put_var(ncid,id,xcord(IRANGE))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xic -")

         status = nf90_inq_varid(ncid,'etac',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","etac_id -")
         do j=1,jmax
            ycord(j) = joff+j
         end do
         status = nf90_put_var(ncid,id,ycord(JRANGE))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","etac -")


         status = nf90_inq_varid(ncid,'xx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xx_id -")
         edges(1) = edges(1) + 1 
         edges(2) = edges(2) + 1
         status = nf90_put_var(ncid,id,xx(-1+IRANGE,-1+JRANGE),start,edges)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xx -")

         status = nf90_inq_varid(ncid,'yx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yx_id -")
         status = nf90_put_var(ncid,id,yx(-1+IRANGE,-1+JRANGE),start,edges)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","yx -")
         edges(1) = edges(1) - 1
         edges(2) = edges(2) - 1

#if 0
         status = nf90_inq_varid(ncid,'lonc',lonc_id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lonc_id -")
         status = nf90_put_vara_double(ncid,latc_id,start,edges,          &
                                     latc(1:imax,1:jmax))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","latc -")

         status = nf90_inq_varid(ncid,'latc',latc_id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","latc_id -")
         status = nf90_put_vara_double(ncid,lonc_id,start,edges,          &
                                     lonc(1:imax,1:jmax))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","lonc -")


         status = nf90_inq_varid(ncid,'convx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                    "save_grid_ncdf()","convx_id -")
         status = nf90_put_vara_double(ncid,convc_id,start,edges,         &
                                        convc(1:imax,1:jmax))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","convc -")
#endif
#endif
      case default

   end select

   if ( save_metrics ) then
      select case (grid_type)

         case (1)

            if ( have_lonlat ) then
               status = nf90_inq_varid(ncid,'latu',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latu -")
               status = nf90_put_var(ncid,id,latu(IRANGE,JRANGE),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latu -")

               status = nf90_inq_varid(ncid,'latv',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latv -")
               status = nf90_put_var(ncid,id,latv(IRANGE,JRANGE),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latv -")

               status = nf90_inq_varid(ncid,'convc',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","convc -")
               status = nf90_put_var(ncid,id,convc(IRANGE,JRANGE),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","convc -")

            end if

         case default

            status = nf90_inq_varid(ncid,'latu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latu -")
            status = nf90_put_var(ncid,id,latu(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latu -")

            status = nf90_inq_varid(ncid,'latv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latv -")
            status = nf90_put_var(ncid,id,latv(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latv -")

            status = nf90_inq_varid(ncid,'dxc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxc -")
            status = nf90_put_var(ncid,id,dxc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxc -")

            status = nf90_inq_varid(ncid,'dyc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyc -")
            status = nf90_put_var(ncid,id,dyc(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyc -")

            status = nf90_inq_varid(ncid,'dxu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxu -")
            status = nf90_put_var(ncid,id,dxu(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxu -")

            status = nf90_inq_varid(ncid,'dyu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyu -")
            status = nf90_put_var(ncid,id,dyu(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyu -")

            status = nf90_inq_varid(ncid,'dxv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxv -")
            status = nf90_put_var(ncid,id,dxv(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxv -")

            status = nf90_inq_varid(ncid,'dyv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyv -")
            status = nf90_put_var(ncid,id,dyv(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyv -")

            status = nf90_inq_varid(ncid,'dxx',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxx -")
            status = nf90_put_var(ncid,id,dxx(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxx -")

            status = nf90_inq_varid(ncid,'dyx',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyx -")
            status = nf90_put_var(ncid,id,dyx(IRANGE,JRANGE),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyx -")

      end select
   end if

!  vertical levels
   if (save3d) then

      select case (vert_cord)
      case (1)
         zname  = 'sigma'
      case (2)
         zname  = 'z'
      case (3,4,5)
         zname  = 'level'
      case default
      end select

      status = nf90_inq_varid(ncid,zname,id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","z_id -")
      select case (vert_cord)
      case (1,2,3,4,5)
         status = nf90_put_var(ncid,id,ga)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","ga -")
      case default
      end select
   endif

!  save bathymetry
   status = nf90_inq_varid(ncid,'bathymetry',id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","bathymetry_id -")
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if(az(i,j) .gt. 0) then
            ws(i,j) = H(i,j)
         else
            ws(i,j) = -10.
         end if
      end do
   end do
   status = nf90_put_var(ncid,id,ws(IRANGE,JRANGE),start,edges)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","bathymetry -")

   if (save_masks) then
      status = nf90_inq_varid(ncid,'t_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask_id")
      status = nf90_put_var(ncid,id,az(IRANGE,JRANGE),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask")

      status = nf90_inq_varid(ncid,'u_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf90_put_var(ncid,id,au(IRANGE,JRANGE),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","u_mask")

      status = nf90_inq_varid(ncid,'v_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf90_put_var(ncid,id,av(IRANGE,JRANGE),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask")
   end if

   return
   end subroutine save_grid_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
