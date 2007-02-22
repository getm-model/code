!$Id: init_grid_ncdf.F90,v 1.3 2007-02-22 08:48:12 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialse grid related variables
!
! !INTERFACE:
   subroutine init_grid_ncdf(ncid,init3d,x_dim,y_dim,z_dim)
!
! !DESCRIPTION:
! This routine creates netCDF variables in an already existing netCDF file
! in define mode with netCDF file-id "{\tt ncid}". All variables are
! related the numerical grid and the bathymetery. If the logical flag
! "{\tt init3d}" evaluates false, no information about the vertical grid
! is initalised (e.g.\ if results from a horizontally integrated run are stored).
! Output arguments are the dimension id's for the netCDF dimensions, which
! may be needed for creating other, not grid related, netCDF variables.
!
! !USES:
   use exceptions
   use output, only: save_masks
   use ncdf_common
   use grid_ncdf
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: xy_exists,latlon_exists
   use domain, only: grid_type,proj_type,vert_cord
   use domain, only: proj_exists
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)            :: ncid
   logical, intent(in)            :: init3d
!
! !INPUT PARAMETERS:
   integer, intent(out)           :: x_dim
   integer, intent(out)           :: y_dim
   integer, intent(out), optional :: z_dim
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: init_grid_ncdf.F90,v $
!  Revision 1.3  2007-02-22 08:48:12  kbk
!  possible to save masks (az, au, av)
!
!  Revision 1.2  2005/04/29 12:45:41  kbk
!  stricter COARDS conforming
!
!  Revision 1.1  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: status
   integer                   :: scalar(1),axisdim(1)
   integer                   :: f2_dims(2),f3_dims(3)
   integer                   :: z_id
   integer                   :: grid_type_id,proj_type_id,vert_cord_id
   integer                   :: proj_lat_id,proj_lon_id,proj_rot_id
   integer                   :: rearth_id
   integer                   :: ioff_id,joff_id
   integer                   :: x0_id,y0_id,lon0_id,lat0_id
   integer                   :: dx_id,dy_id,dlon_id,dlat_id
   integer                   :: xc_id,yc_id
   integer                   :: lonc_id,latc_id
   integer                   :: convc_id
   integer                   :: t_mask_id,u_mask_id,v_mask_id
   integer                   :: bathymetry_id
   REALTYPE                  :: fv,mv,vr(2)
   character(32)             :: xname,yname,zname
   character(32)             :: xunits,yunits,zunits
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"

!  check if function is called with correct arguments
   if (present(z_dim) .and. ( .not. init3d)) then
      call getm_error("init_grid_ncdf()",                               &
                      "Dummy argument 'z_dim' illegally passed.")
   endif

   if (( .not. present(z_dim)) .and. init3d) then
      call getm_error("init_grid_ncdf()",                               &
                      "Dummy argument 'z_dim' missing.")
   endif

!  length of netCDF dimensions
   xlen = imax-imin+1
   ylen = jmax-jmin+1
   zlen = kmax+1


!  some grid-specific settings
   select case (grid_type)
      case (1)
         xname    = 'xc'
         yname    = 'yc'
         xunits   = 'm'
         yunits   = 'm'
      case (2)
         xname    = 'lonc'
         yname    = 'latc'
         xunits   = 'degrees_east'
         yunits   = 'degrees_north'
      case (3)
         xname    = 'xic'
         yname    = 'etac'
         xunits   = ' '
         yunits   = ' '
      case (4)
         xname    = 'xic'
         yname    = 'etac'
         xunits   = ' '
         yunits   = ' '
      case default
   end select

   if (init3d) then
      select case (vert_cord)
      case (1)
         zname  = 'sigma'
         zunits = 'sigma_level'
      case (2)
         zname  = 'z'
         zunits = 'm'
      case (3)
         zname  = 'level'
         zunits = 'level'
      case default
      end select
   endif


!  create dimensions
   status = nf_def_dim(ncid,xname,xlen,x_dim)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","x_dim -")

   status = nf_def_dim(ncid,yname,ylen,y_dim)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","y_dim -")

   if (init3d) then
      status = nf_def_dim(ncid,zname,zlen,z_dim)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","z_dim -")
   endif


!  netCDF dimension vectors
   f2_dims(2)= y_dim
   f2_dims(1)= x_dim


   if (init3d) then
      f3_dims(3)= z_dim
      f3_dims(2)= y_dim
      f3_dims(1)= x_dim
   endif


!  horizontal grid types
   scalar(1) = 0
   status = nf_def_var(ncid,'grid_type',NF_INT,0,scalar,grid_type_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","grid_type -")

   status = nf_put_att_text(ncid,grid_type_id,'option_1_',              &
                         len_trim('cartesian'),'cartesian')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","cartesian -")

   status = nf_put_att_text(ncid,grid_type_id,'option_2_',              &
                         len_trim('spherical'),'spherical')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","spherical -")

   status = nf_put_att_text(ncid,grid_type_id,'option_3_',              &
                         len_trim('curvilinear'),'curvilinear')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","curvilinear -")

   status = nf_put_att_text(ncid,grid_type_id,'option_4_',              &
                         len_trim('spherilinear'), 'spherilinear')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","spherilinear -")


!  vertical grid types
   if (init3d) then
      status = nf_def_var(ncid,'vert_cord',NF_INT,0,scalar,vert_cord_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","vert_cord -")

      status = nf_put_att_text(ncid,vert_cord_id,'option_1_',           &
                               len_trim('sigma'),'sigma')
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","sigma -")

      status = nf_put_att_text(ncid,vert_cord_id,'option_2_',           &
                               len_trim('z'),'z')
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","z -")

      status = nf_put_att_text(ncid,vert_cord_id,'option_3_',           &
                              len_trim('s'),'s')
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","s -")
   endif

!  geographic projection types
   status = nf_def_var(ncid,'proj_type',NF_INT,0,scalar,proj_type_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","proj_type -")


   status = nf_put_att_text(ncid,proj_type_id,'option_1_',              &
                         len_trim('mercator'),'mercator')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","mercator -")

   status = nf_put_att_text(ncid,proj_type_id,'option_2_',              &
                         len_trim('stereographic'),'stereographic')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","stereographic -")

   status = nf_put_att_text(ncid,proj_type_id,'option_3_',              &
                         len_trim('lambert'),'lambert')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","lambert -")

   status = nf_put_att_text(ncid,proj_type_id,'option_99_',             &
                         len_trim('unknown'),'unknown')
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","unkown -")

!  geographic projection variables
   status = nf_def_var(ncid,'proj_lat',NF_REAL,0,scalar,proj_lat_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","proj_lat -")
   call set_attributes(ncid,proj_lat_id,long_name='reference latitude')

   status = nf_def_var(ncid,'proj_lon',NF_REAL,0,scalar,proj_lon_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","proj_lon -")
   call set_attributes(ncid,proj_lon_id,long_name='reference longitude')

   status = nf_def_var(ncid,'proj_rot',NF_REAL,0,scalar,proj_rot_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","proj_rot -")
   call set_attributes(ncid,proj_rot_id,long_name='map rotation')

   status = nf_def_var(ncid,'rearth',NF_REAL,0,scalar,rearth_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","rearth -")
   call set_attributes(ncid,rearth_id,long_name='radius of earth')



!  offset for parallel runs
   status = nf_def_var(ncid,'ioff',NF_INT,0,scalar,ioff_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","ioff -")
   call set_attributes(ncid,ioff_id,long_name='index offset (i)')

   status = nf_def_var(ncid,'joff',NF_INT,0,scalar,joff_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","joff -")
   call set_attributes(ncid,joff_id,long_name='index offset (j)')


!  grid-specific settings
   select case (grid_type)
      case (1)

!        grid spacing
         status = nf_def_var(ncid,'dx',NF_REAL,0,scalar,dx_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dx -")
         call set_attributes(ncid,dx_id,units='m',                      &
                             long_name='grid spacing (x)')

         status = nf_def_var(ncid,'dy',NF_REAL,0,scalar,dy_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dy -")
         call set_attributes(ncid,dy_id,units='m',                      &
                             long_name='grid spacing (y)')

!        global offset
         status = nf_def_var(ncid,'x0',NF_REAL,0,scalar,x0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "init_grid_ncdf()","x0 -")
         call set_attributes(ncid,x0_id,units='m',                      &
                             long_name='offset (x)')

         status = nf_def_var(ncid,'y0',NF_REAL,0,scalar,y0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","y0 -")
         call set_attributes(ncid,y0_id,units='m',                      &
                             long_name='offset (y)')

!        coordinate variables
         axisdim(1) = x_dim
         status = nf_def_var(ncid,'xc',NF_REAL,1,axisdim,xc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                       "init_grid_ncdf()","xc -")
         call set_attributes(ncid,xc_id,units=trim(xunits))

         axisdim(1) = y_dim
         status = nf_def_var(ncid,'yc',NF_REAL,1,axisdim,yc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","yc -")
         call set_attributes(ncid,yc_id,units=trim(yunits))

      case (2)

!        grid spacing
         status = nf_def_var(ncid,'dlon',NF_REAL,0,scalar,dlon_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dlon -")
         call set_attributes(ncid,dlon_id,units=trim(xunits),           &
                             long_name='grid spacing (lon)')

         status = nf_def_var(ncid,'dlat',NF_REAL,0,scalar,dlat_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dlat -")
         call set_attributes(ncid,dlat_id,units=trim(yunits),           &
                             long_name='grid spacing (lat)')

!        global offset
         status = nf_def_var(ncid,'lon0',NF_REAL,0,scalar,lon0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","lon0 -")
         call set_attributes(ncid,lon0_id,units=trim(xunits),           &
                             long_name='offset (lon)')

         status = nf_def_var(ncid,'lat0',NF_REAL,0,scalar,lat0_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","lat0 -")
         call set_attributes(ncid,lat0_id,units=trim(yunits),           &
                             long_name='offset (lat)')

!        coordinate variables
         axisdim(1) = x_dim
         status = nf_def_var(ncid,'lonc',NF_REAL,1,axisdim,lonc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","lonc -")
         call set_attributes(ncid,lonc_id,units=trim(xunits))

         axisdim(1) = y_dim
         status = nf_def_var(ncid,'latc',NF_REAL,1,axisdim,latc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","latc -")
         call set_attributes(ncid,latc_id,units=trim(yunits))
      case (3)
!        no netcdf support for coordinate variables
!        for irregularly spaced coordinates
      case (4)
!        no netcdf support for coordinate variables
!        for irregularly spaced coordinates
      case default
   end select

   if (init3d) then
      axisdim(1) = z_dim
      status = nf_def_var(ncid,zname,NF_REAL,1,axisdim,z_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                    "init_grid_ncdf()","z -")
      call set_attributes(ncid,z_id,units=zunits)
   endif



!  x and y positions of T-points
!  for non-Cartesian grids
   if (grid_type .ne. 1) then

      if (xy_exists) then

         status = nf_def_var(ncid,'xc',NF_REAL,2,f2_dims,xc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","xc -")
         fv = xy_missing
         mv = xy_missing
         vr(1) = -1.e8
         vr(2) =  1.e8
         call set_attributes(ncid,xc_id,  &
                            long_name='x-position',units='m',          &
                            FillValue=fv,missing_value=mv,valid_range=vr)

         status = nf_def_var(ncid,'yc',NF_REAL,2,f2_dims,yc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","yc -")
         fv = xy_missing
         mv = xy_missing
         vr(1) = -1.e8
         vr(2) =  1.e8
         call set_attributes(ncid,yc_id,  &
                            long_name='y-position',units='m',           &
                            FillValue=fv,missing_value=mv,valid_range=vr)

      endif

   endif

!  bathymetry
   status = nf_def_var(ncid,'bathymetry',NF_REAL,2,f2_dims,bathymetry_id)
   if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","bathymetry -")
   fv = h_missing
   mv = h_missing
   vr(1) = -5.
   vr(2) = 4000.
   call set_attributes(ncid,bathymetry_id,                              &
                       long_name='bathymetry',units='m',                &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   if (save_masks) then
      status = nf_def_var(ncid,'t_mask',NF_INT,2,f2_dims,t_mask_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","t_mask")
      call set_attributes(ncid,t_mask_id,long_name='mask for T-points')

      status = nf_def_var(ncid,'u_mask',NF_INT,2,f2_dims,u_mask_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","u_mask")
      call set_attributes(ncid,u_mask_id,long_name='mask for U-points')

      status = nf_def_var(ncid,'v_mask',NF_INT,2,f2_dims,v_mask_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","v_mask")
      call set_attributes(ncid,v_mask_id,long_name='mask for V-points')
   end if

!  lat,lon positions of T-points and grid rotation
!  for non-spherical grids
   if (grid_type .ne. 2) then

      if (latlon_exists) then

         status = nf_def_var(ncid,'latc',NF_REAL,2,f2_dims,latc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","latc -")
         fv = latlon_missing
         mv = latlon_missing
         vr(1) = -90.
         vr(2) =  90.
         call set_attributes(ncid,latc_id,  &
                            long_name='latitude',units='degrees_north', &
                            FillValue=fv,missing_value=mv,valid_range=vr)

         status = nf_def_var(ncid,'lonc',NF_REAL,2,f2_dims,lonc_id)
         if (status .ne. NF_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","lonc -")
         fv = latlon_missing
         mv = latlon_missing
         vr(1) = -180.
         vr(2) =  180.
         call set_attributes(ncid,lonc_id,  &
                            long_name='longitude',units='degrees_east', &
                            FillValue=fv,missing_value=mv,valid_range=vr)

      endif

!     angle of rotation between local grid axes and N-S axes
      status = nf_def_var(ncid,'convc',NF_REAL,2,f2_dims,convc_id)
      if (status .ne. NF_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","convc -")
      fv = conv_missing
      mv = conv_missing
      vr(1) = -180.
      vr(2) =  180.
      call set_attributes(ncid,convc_id,  &
                         long_name='grid rotation',units='degrees',     &
                         FillValue=fv,missing_value=mv,valid_range=vr)
   endif

   return
   end subroutine init_grid_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
