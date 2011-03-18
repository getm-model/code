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
   use domain, only: xxcord,yxcord
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
! set _WRITE_HALOS_ via Makefile
#ifdef _WRITE_HALOS_
   LEVEL3 'include HALOs in NetCDF output files'
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
         status = nf90_put_var(ncid,id,xcord(_IRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xc -")

         status = nf90_inq_varid(ncid,'yc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yc_id -")
         status = nf90_put_var(ncid,id,ycord(_JRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","yc -")

         if ( have_lonlat ) then

            status = nf90_inq_varid(ncid,'lonc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","lonc -")
            status = nf90_put_var(ncid,id,lonc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","lonc -")

            status = nf90_inq_varid(ncid,'latc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latc -")
            status = nf90_put_var(ncid,id,latc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latc -")

            status = nf90_inq_varid(ncid,'convc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","convc -")
            status = nf90_put_var(ncid,id,convc(_2D_W_),start,edges)
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
         status = nf90_put_var(ncid,id,xcord(_IRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","lonc -")

         status = nf90_inq_varid(ncid,'latc',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","latc_id -")
         status = nf90_put_var(ncid,id,ycord(_JRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","latc -")
         if ( have_xy ) then

            status = nf90_inq_varid(ncid,'xc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","xc -")
            status = nf90_put_var(ncid,id,xc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","xc -")

            status = nf90_inq_varid(ncid,'yc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","yc -")
            status = nf90_put_var(ncid,id,yc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","yc -")

         end if

      case (3,4)

!        pseudo coordinate variables - T-points
         status = nf90_inq_varid(ncid,'xic',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xic_id -")
         status = nf90_put_var(ncid,id,xcord(_IRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xic -")

         status = nf90_inq_varid(ncid,'etac',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","etac_id -")
         status = nf90_put_var(ncid,id,ycord(_JRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","etac -")

!        pseudo coordinate variables - X-points
         status = nf90_inq_varid(ncid,'xix',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xix_id -")
         status = nf90_put_var(ncid,id,xxcord(-1+_IRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xix -")

         status = nf90_inq_varid(ncid,'etax',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","etax_id -")
         status = nf90_put_var(ncid,id,yxcord(-1+_JRANGE_))
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","etax -")

!        positions of vertices
         status = nf90_inq_varid(ncid,'xx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xx_id -")
         edges(1) = edges(1) + 1 
         edges(2) = edges(2) + 1
         status = nf90_put_var(ncid,id,xx(-1+_IRANGE_,-1+_JRANGE_),start,edges)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xx -")

         status = nf90_inq_varid(ncid,'yx',id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yx_id -")
         status = nf90_put_var(ncid,id,yx(-1+_IRANGE_,-1+_JRANGE_),start,edges)
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
      case default

   end select

   if ( save_metrics ) then
      select case (grid_type)

         case (1)

            if ( have_lonlat ) then
               status = nf90_inq_varid(ncid,'latu',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latu -")
               status = nf90_put_var(ncid,id,latu(_2D_W_),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latu -")

               status = nf90_inq_varid(ncid,'latv',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latv -")
               status = nf90_put_var(ncid,id,latv(_2D_W_),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","latv -")

               status = nf90_inq_varid(ncid,'convc',id)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","convc -")
               status = nf90_put_var(ncid,id,convc(_2D_W_),start,edges)
               if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                             "save_grid_ncdf()","convc -")

            end if

         case default

            status = nf90_inq_varid(ncid,'latu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latu -")
            status = nf90_put_var(ncid,id,latu(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latu -")

            status = nf90_inq_varid(ncid,'latv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latv -")
            status = nf90_put_var(ncid,id,latv(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","latv -")

            status = nf90_inq_varid(ncid,'dxc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxc -")
            status = nf90_put_var(ncid,id,dxc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxc -")

            status = nf90_inq_varid(ncid,'dyc',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyc -")
            status = nf90_put_var(ncid,id,dyc(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyc -")

            status = nf90_inq_varid(ncid,'dxu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxu -")
            status = nf90_put_var(ncid,id,dxu(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxu -")

            status = nf90_inq_varid(ncid,'dyu',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyu -")
            status = nf90_put_var(ncid,id,dyu(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyu -")

            status = nf90_inq_varid(ncid,'dxv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxv -")
            status = nf90_put_var(ncid,id,dxv(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxv -")

            status = nf90_inq_varid(ncid,'dyv',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyv -")
            status = nf90_put_var(ncid,id,dyv(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyv -")

            status = nf90_inq_varid(ncid,'dxx',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxx -")
            status = nf90_put_var(ncid,id,dxx(_2D_W_),start,edges)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dxx -")

            status = nf90_inq_varid(ncid,'dyx',id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                          "save_grid_ncdf()","dyx -")
            status = nf90_put_var(ncid,id,dyx(_2D_W_),start,edges)
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
   status = nf90_put_var(ncid,id,ws(_2D_W_),start,edges)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","bathymetry -")

   if (save_masks) then
      status = nf90_inq_varid(ncid,'t_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask_id")
      status = nf90_put_var(ncid,id,az(_2D_W_),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask")

      status = nf90_inq_varid(ncid,'u_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf90_put_var(ncid,id,au(_2D_W_),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","u_mask")

      status = nf90_inq_varid(ncid,'v_mask',id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf90_put_var(ncid,id,av(_2D_W_),start,edges)
      if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask")
   end if

   status = nf90_sync(ncid)
   if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                     "save_grid_ncdf()","syncing")

   return
   end subroutine save_grid_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
