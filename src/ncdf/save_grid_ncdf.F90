!$Id: save_grid_ncdf.F90,v 1.6 2007-10-16 07:14:35 kbk Exp $
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
   use output, only: save_masks
   use grid_ncdf
   use domain, only: imin,imax,jmin,jmax
   use domain, only: ioff,joff
   use domain, only: x0,y0,dx,dy
   use domain, only: lon0,lat0,dlon,dlat
   use domain, only: xy_exists,xc,yc
   use domain, only: latlon_exists,latc,lonc,convc
   use domain, only: grid_type,proj_type,vert_cord
   use domain, only: proj_exists
   use domain, only: proj_lon,proj_lat,proj_rot
   use domain, only: rearth
   use domain, only: h,ga
   use domain, only: az,au,av

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

   integer                   :: z_id

   integer                   :: grid_type_id,vert_cord_id
   integer                   :: proj_type_id
   integer                   :: proj_lat_id,proj_lon_id,proj_rot_id
   integer                   :: rearth_id
   integer                   :: ioff_id,joff_id
   integer                   :: x0_id,y0_id,lon0_id,lat0_id
   integer                   :: dx_id,dy_id,dlon_id,dlat_id
   integer                   :: xc_id,yc_id
   integer                   :: lonc_id,latc_id
   integer                   :: xic_id,etac_id
   integer                   :: convc_id
   integer                   :: bathymetry_id,id

   character(32)             :: zname


   REAL_4B, dimension(:), allocatable :: ws
   REALTYPE, dimension(:), allocatable :: cord
!
!EOP
!------------------------------------------------------------------------
!BOC
   include "netcdf.inc"

!  allocate work space
   allocate(ws(xlen*ylen),stat=status)
   if (status .ne. 0) call getm_error("save_grid_ncdf()",               &
                                      "error allocating ws")


!----------- get netCDF ID's --------------------------------------------

!  get grid and projection id's
   status = nf_inq_varid(ncid,'grid_type',grid_type_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","grid_type_id -")

   if (save3d) then
      status = nf_inq_varid(ncid,'vert_cord',vert_cord_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","vert_cord_id -")
   endif

   status = nf_inq_varid(ncid,'proj_type',proj_type_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_type_id -")

   status = nf_inq_varid(ncid,'proj_lat',proj_lat_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_lat_id -")

   status = nf_inq_varid(ncid,'proj_lon',proj_lon_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_lon_id -")

   status = nf_inq_varid(ncid,'proj_rot',proj_rot_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_rot_id -")

   status = nf_inq_varid(ncid,'rearth',rearth_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","rearth_id -")


! offset for parallel runs
   status = nf_inq_varid(ncid,'ioff',ioff_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","ioff_id -")

   status = nf_inq_varid(ncid,'joff',joff_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","joff_id -")


!  grid-specific settings
   select case (grid_type)
      case (1)

!        grid spacing
         status = nf_inq_varid(ncid,'dx',dx_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dx_id -")

         status = nf_inq_varid(ncid,'dy',dy_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dy_id -")


!        global set-off
         status = nf_inq_varid(ncid,'x0',x0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","x0_id -")

         status = nf_inq_varid(ncid,'y0',y0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","y0_id -")

!        coordinate variables
         status = nf_inq_varid(ncid,'xc',xc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xc_id -")

         status = nf_inq_varid(ncid,'yc',yc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yc_id -")

      case (2)

!        grid spacing
         status = nf_inq_varid(ncid,'dlon',dlon_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dlon_id -")

         status = nf_inq_varid(ncid,'dlat',dlat_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","dlat_id -")

!        global set-off
         status = nf_inq_varid(ncid,'lon0',lon0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lon0_id -")

         status = nf_inq_varid(ncid,'lat0',lat0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lat0_id -")

!        coordinate variables
         status = nf_inq_varid(ncid,'lonc',lonc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lonc_id -")

         status = nf_inq_varid(ncid,'latc',latc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","latc_id -")

      case (3,4)
!        pseudo coordinate variables
         status = nf_inq_varid(ncid,'xic',xic_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xic_id -")

         status = nf_inq_varid(ncid,'etac',etac_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","etac_id -")
      case default
   end select

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

      status = nf_inq_varid(ncid,zname,z_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","z_id -")
   endif

   if (grid_type.ne.1) then
      if (xy_exists) then
         status = nf_inq_varid(ncid,'xc',xc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","xc_id -")

         status = nf_inq_varid(ncid,'yc',yc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","yc_id -")
      endif
   endif

   if (grid_type.ne.2) then
      if (latlon_exists) then
         status = nf_inq_varid(ncid,'lonc',lonc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","lonc_id -")

         status = nf_inq_varid(ncid,'latc',latc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","latc_id -")

      endif

      status = nf_inq_varid(ncid,'convc',convc_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                    "save_grid_ncdf()","convc_id -")
   endif


   status = nf_inq_varid(ncid,'bathymetry',bathymetry_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","bathymetry_id -")




!----------- save netCDF variables  -------------------------------------


!  grid types
   status = nf_put_var_int(ncid,grid_type_id,grid_type)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","grid_type -")

   if (save3d) then
      status = nf_put_var_int(ncid,vert_cord_id,vert_cord)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                    "save_grid_ncdf()","vert_cord -")
   endif

! projection types
   status = nf_put_var_int(ncid,proj_type_id,proj_type)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_type -")

   status = nf_put_var_double(ncid,proj_lat_id,proj_lat)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_lat -")

   status = nf_put_var_double(ncid,proj_lon_id,proj_lon)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_lon -")

   status = nf_put_var_double(ncid,proj_rot_id,proj_rot)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","proj_rot -")

   status = nf_put_var_double(ncid,rearth_id,rearth)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","rearth -")

!  domain offset for parallel runs
   status = nf_put_var_int(ncid,ioff_id,ioff)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","ioff -")

   status = nf_put_var_int(ncid,joff_id,joff)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","joff -")

!  save coordinate information
!  Using the F77 NetCDF interface it is necessary to make a copy of the
!  coordinate variables to a 1-based vector. The cord variable is 
!  allocated here and used further down.
   allocate(cord(max(imax,jmax)),stat=status)
   if (status .ne. 0) call getm_error("save_grid_ncdf()",               &
                                      "error allocating cord")

   select case (grid_type)
   case (1)

!     global offset
      status = nf_put_var_double(ncid,dx_id,dx)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","dx -")

      status = nf_put_var_double(ncid,dy_id,dy)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","dy -")

!     grid spacing
      status = nf_put_var_double(ncid,x0_id,x0)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","x0 -")

      status = nf_put_var_double(ncid,y0_id,y0)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","y0 -")

!     coordinate variables
#if 1
      cord=xc(1:imax,1)
      status = nf_put_var_double(ncid,xc_id,cord)
#else
      status = nf_put_var_double(ncid,xc_id,xc(1:imax,1))
#endif
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","xc -")
#if 1
      cord=yc(1,1:jmax)
      status = nf_put_var_double(ncid,yc_id,cord)
#else
      status = nf_put_var_double(ncid,yc_id,yc(1,1:jmax))
#endif
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","yc -")
   case (2)

!     global offset
      status = nf_put_var_double(ncid,dlat_id,dlat)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","dlat -")

      status = nf_put_var_double(ncid,dlon_id,dlon)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","dlon -")

!     grid spacing
      status = nf_put_var_double(ncid,lat0_id,lat0)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","lat0 -")

      status = nf_put_var_double(ncid,lon0_id,lon0)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","lon0 -")

!     coordinate variables
#if 1
      cord=lonc(1:imax,1)
      status = nf_put_var_double(ncid,lonc_id,cord)
#else
      status = nf_put_var_double(ncid,lonc_id,lonc(1:imax,1))
#endif
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","lonc -")
#if 1
      cord=latc(1,1:jmax)
      status = nf_put_var_double(ncid,latc_id,cord)
#else
      status = nf_put_var_double(ncid,latc_id,latc(1,1:jmax))
#endif
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","latc -")
   case (3,4)
!     pseudo coordinate variables
      do i=1,imax
         cord(i) = ioff+i
      end do
      status = nf_put_var_double(ncid,xic_id,cord)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","xic -")
      do j=1,jmax
         cord(j) = joff+j
      end do
      status = nf_put_var_double(ncid,etac_id,cord)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","etac -")
   case default
   end select


!  vertical levels
   if (save3d) then
      select case (vert_cord)
      case (1,2,3,4,5)
         status = nf_put_var_double(ncid,z_id,ga)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "save_grid_ncdf()","ga -")
      case default
      end select
   endif

!  Set netCDF slice information
   start(1) = 1
   start(2) = 1
   edges(1) = xlen
   edges(2) = ylen

!  save bathymetry
   call cnv_2d(imin,jmin,imax,jmax,az,H,h_missing,                      &
               imin,jmin,imax,jmax,ws)
   status = nf_put_vara_real(ncid,bathymetry_id,start,edges,ws)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "save_grid_ncdf()","bathymetry -")

   if (save_masks) then
      status = nf_inq_varid(ncid,'t_mask',id)
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask_id")
      status = nf_put_vara_int(ncid,id,start,edges,az(imin:imax,jmin:jmax))
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","t_mask")

      status = nf_inq_varid(ncid,'u_mask',id)
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf_put_vara_int(ncid,id,start,edges,au(imin:imax,jmin:jmax))
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","u_mask")

      status = nf_inq_varid(ncid,'v_mask',id)
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask_id")
      status = nf_put_vara_int(ncid,id,start,edges,av(imin:imax,jmin:jmax))
      if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","v_mask")
   end if

!  Save x and y position of T-points
   if (grid_type .ne. 1) then

      if (xy_exists) then

         status = nf_put_vara_double(ncid,xc_id,start,edges,            &
                                     xc(1:imax,1:jmax))
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","xc -")

         status = nf_put_vara_double(ncid,yc_id,start,edges,            &
                                     yc(1:imax,1:jmax))
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","yc -")

      endif

   endif

!  Save lat and lon of T-points
   if (grid_type .ne. 2) then

      if (latlon_exists) then

         status = nf_put_vara_double(ncid,latc_id,start,edges,          &
                                     latc(1:imax,1:jmax))
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","latc -")

         status = nf_put_vara_double(ncid,lonc_id,start,edges,          &
                                     lonc(1:imax,1:jmax))
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "save_grid_ncdf()","lonc -")

      endif

!     save grid rotation
      status = nf_put_vara_double(ncid,convc_id,start,edges,            &
                                     convc(1:imax,1:jmax))
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "save_grid_ncdf()","convc -")

   endif

!  reclaim storage
   if(allocated(ws)) deallocate(ws,stat=status)
   if (status .ne. 0) call getm_error("save_grid_ncdf()",               &
                                      "error deallocating ws")

   return
   end subroutine save_grid_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
