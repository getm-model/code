!$Id: init_grid_ncdf.F90,v 1.7 2009-09-23 10:11:48 kb Exp $
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
   use netcdf
   use ncdf_common
   use grid_ncdf
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: grid_type,vert_cord
   use domain, only: have_lonlat,have_xy
   use output, only: save_metrics,save_masks
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
!  Revision 1.7  2009-09-23 10:11:48  kb
!  rewrite of grid-initialisation, optional grid info saved to file, -DSAVE_HALO, updated documentation
!
!  Revision 1.6  2009-03-13 14:44:14  kb
!  grid information in NF_DOUBLE
!
!  Revision 1.5  2007-10-16 07:14:35  kbk
!  pseudo coordinate variables for curvi-linear grids
!
!  Revision 1.4  2007-03-30 13:11:00  hb
!  Use of adaptive and hybrid vertical coordinates technically enabled
!
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
   integer                   :: id
   integer                   :: axisdim(1)
   integer                   :: f2_dims(2)
   REALTYPE                  :: fv,mv,vr(2)
   character(32)             :: xname,yname,zname
   character(32)             :: xunits,yunits,zunits
!EOP
!-----------------------------------------------------------------------
!BOC

!  check if function is called with correct arguments
   if (present(z_dim) .and. ( .not. init3d)) then
      call getm_error("init_grid_ncdf()",                               &
                      "Dummy argument 'z_dim' illegally passed.")
   endif

   if (( .not. present(z_dim)) .and. init3d) then
      call getm_error("init_grid_ncdf()",                               &
                      "Dummy argument 'z_dim' missing.")
   endif

! length of netCDF dimensions
! set INCLUDE_HALOS via Makefile
#ifdef INCLUDE_HALOS
   xlen = (imax+HALO)-(imin-HALO)+1
   ylen = (jmax+HALO)-(jmin-HALO)+1
   zlen = kmax+1
#else
   xlen = imax-imin+1
   ylen = jmax-jmin+1
   zlen = kmax+1
#endif

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
      case (3,4,5)
         zname  = 'level'
         zunits = 'level'
      case default
      end select
   endif

!  create dimensions

   status = nf90_def_dim(ncid,xname,xlen,xc_dim)
   x_dim=xc_dim
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","xc_dim -")

   status = nf90_def_dim(ncid,yname,ylen,yc_dim)
   y_dim=yc_dim
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","yc_dim -")
   select case (vert_cord)
      case(3,4)
         status = nf90_def_dim(ncid,'kurt',xlen+1,xx_dim)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","xx_dim -")

         status = nf90_def_dim(ncid,'egon',ylen+1,yx_dim)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","yx_dim -")
      case default
   end select

   if (init3d) then
      status = nf90_def_dim(ncid,zname,zlen,z_dim)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","z_dim -")
   endif

!  netCDF dimension vectors
   f2_dims(1)= x_dim
   f2_dims(2)= y_dim

!  horizontal grid types
   status = nf90_def_var(ncid,'grid_type',NF90_INT,id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","grid_type")

!  vertical grid types
   if (init3d) then
      status = nf90_def_var(ncid,'vert_cord',NF90_INT,id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","vert_cord")

!KB      status = nf90_put_att_text(ncid,id,'option_1_',           &
!KB                               len_trim('sigma'),'sigma')
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                     "init_grid_ncdf()","sigma -")
   endif

!  offset for parallel runs
   status = nf90_def_var(ncid,'ioff',NF90_INT,id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","ioff -")
   call set_attributes(ncid,id,long_name='index offset (i)')

   status = nf90_def_var(ncid,'joff',NF90_INT,id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","joff -")
   call set_attributes(ncid,id,long_name='index offset (j)')


!  grid-specific settings
   select case (grid_type)
      case (1)

!        grid spacing
         status = nf90_def_var(ncid,'dx',NF90_DOUBLE,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dx -")
         call set_attributes(ncid,id,units='m',                         &
                             long_name='grid spacing (x)')

         status = nf90_def_var(ncid,'dy',NF90_DOUBLE,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dy -")
         call set_attributes(ncid,id,units='m',                      &
                             long_name='grid spacing (y)')

!        coordinate variables
         status = nf90_def_var(ncid,'xc',NF90_DOUBLE,(/ x_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                       "init_grid_ncdf()","xc -")
         call set_attributes(ncid,id,units=trim(xunits))

         status = nf90_def_var(ncid,'yc',NF90_DOUBLE,(/ y_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","yc -")
         call set_attributes(ncid,id,units=trim(yunits))

!        longtitude, latitude information if present
         if ( have_lonlat ) then

!           lonc
            status = nf90_def_var(ncid,'lonc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                           "init_grid_ncdf()","lonc -")
            fv = latlon_missing
            mv = latlon_missing
            vr(1) = -180.
            vr(2) =  180.
            call set_attributes(ncid,id,  &
                               long_name='longitude',units='degrees_east', &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           latc
            status = nf90_def_var(ncid,'latc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                           "init_grid_ncdf()","latc -")
            fv = latlon_missing
            mv = latlon_missing
            vr(1) = -90.
            vr(2) =  90.
            call set_attributes(ncid,id,  &
                               long_name='latitude',units='degrees_north', &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           angle of rotation between local grid axes and N-S axes
            status = nf90_def_var(ncid,'convc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","convc -")
            fv = conv_missing
            mv = conv_missing
            vr(1) = -180.
            vr(2) =  180.
            call set_attributes(ncid,id,  &
                               long_name='grid rotation',units='degrees',     &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           latitude of U-points
            status = nf90_def_var(ncid,'latu',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","latu -")
            fv = conv_missing
            mv = conv_missing
            vr(1) = -90.
            vr(2) =  90.
            call set_attributes(ncid,id,  &
                               long_name='latu',units='degrees',           &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           latitude of V-points
            status = nf90_def_var(ncid,'latv',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","latv -")
            fv = conv_missing
            mv = conv_missing
            vr(1) = -90.
            vr(2) =  90.
            call set_attributes(ncid,id,  &
                               long_name='latv',units='degrees',           &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)
         end if

      case (2)

!        grid spacing
         status = nf90_def_var(ncid,'dlon',NF90_DOUBLE,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dlon -")
         call set_attributes(ncid,id,units=trim(xunits),           &
                             long_name='grid spacing (lon)')

         status = nf90_def_var(ncid,'dlat',NF90_DOUBLE,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","dlat -")
         call set_attributes(ncid,id,units=trim(yunits),           &
                             long_name='grid spacing (lat)')

!        coordinate variables
         status = nf90_def_var(ncid,'lonc',NF90_DOUBLE,(/ x_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","lonc -")
         call set_attributes(ncid,id,units=trim(xunits))

         status = nf90_def_var(ncid,'latc',NF90_DOUBLE,(/ y_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","latc -")
         call set_attributes(ncid,id,units=trim(yunits))

!        x, y information if present
         if ( have_xy ) then

!           xc
            status = nf90_def_var(ncid,'xc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                           "init_grid_ncdf()","lonc -")
            fv = latlon_missing
            mv = latlon_missing
            vr(1) = -180.
            vr(2) =  180.
            call set_attributes(ncid,id,  &
                               long_name='longitude',units='degrees_east', &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           yc
            status = nf90_def_var(ncid,'yc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","latc -")
            fv = latlon_missing
            mv = latlon_missing
            vr(1) = -90.
            vr(2) =  90.
            call set_attributes(ncid,id,  &
                               long_name='latitude',units='degrees_north', &
                               netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

         end if


#if 1
      case (3,4)
!        pseudo coordinate variables
         axisdim(1) = x_dim
         status = nf90_def_var(ncid,'xic',NF90_DOUBLE,axisdim,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","xic -")

         axisdim(1) = y_dim
         status = nf90_def_var(ncid,'etac',NF90_DOUBLE,axisdim,id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","etac -")

         status = nf90_def_var(ncid,'xx',NF90_DOUBLE,(/ xx_dim, yx_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","xc -")
         fv = xy_missing
         mv = xy_missing
         vr(1) = -1.e8
         vr(2) =  1.e8
         call set_attributes(ncid,id,&
                            long_name='xx-position',units='m',          &
                            netcdf_real=NF90_DOUBLE, &
                            FillValue=fv,missing_value=mv,valid_range=vr)

         status = nf90_def_var(ncid,'yx',NF90_DOUBLE,(/ xx_dim, yx_dim /),id)
         if (status .ne. NF90_NOERR) call netcdf_error(status,            &
                                        "init_grid_ncdf()","yc -")
         fv = xy_missing
         mv = xy_missing
         vr(1) = -1.e8
         vr(2) =  1.e8
         call set_attributes(ncid,id,  &
                            long_name='yx-position',units='m',           &
                            netcdf_real=NF90_DOUBLE, &
                            FillValue=fv,missing_value=mv,valid_range=vr)

#endif
      case default
   end select

!  metrics information
   if ( save_metrics ) then
      select case (grid_type)
         case (1)

         case default

!           latitude of U-points
            status = nf90_def_var(ncid,'latu',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","latu -")
            fv = conv_missing
            mv = conv_missing
            vr(1) = -90.
            vr(2) =  90.
            call set_attributes(ncid,id,  &
                               long_name='latitude for U-points', &
                               units='degrees',netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           latitude of V-points
            status = nf90_def_var(ncid,'latv',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,          &
                                           "init_grid_ncdf()","latv -")
            call set_attributes(ncid,id,  &
                               long_name='latitude for V-points', &
                               units='degrees', netcdf_real=NF90_DOUBLE, &
                               FillValue=fv,missing_value=mv,valid_range=vr)

!           metric coefficients
            fv = -999.
            mv = -999.
            status = nf90_def_var(ncid,'dxc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dxc -")
            call set_attributes(ncid,id,  &
                               long_name='dx for T-points',units='m',     &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dyc',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dyc -")
            call set_attributes(ncid,id,  &
                               long_name='dy for T-points',units='m',     &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)
 
            status = nf90_def_var(ncid,'dxu',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dxu -")
            call set_attributes(ncid,id,  &
                               long_name='dx for U-points',units='m',     &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dyu',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dyu -")
            call set_attributes(ncid,id,  &
                               long_name='dyu for U-points',units='m',    &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dxv',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dxv -")
            call set_attributes(ncid,id,  &
                               long_name='dx for V-points',units='m',     &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dyv',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dyv -")
            call set_attributes(ncid,id,  &
                               long_name='dy for V-points',units='m',    &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dxx',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dxx -")
            call set_attributes(ncid,id,  &
                               long_name='dx for X-points',units='m',     &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)

            status = nf90_def_var(ncid,'dyx',NF90_DOUBLE,f2_dims,id)
            if (status .ne. NF90_NOERR) call netcdf_error(status,         &
                                           "init_grid_ncdf()","dyx -")
            call set_attributes(ncid,id,  &
                               long_name='dy for X-points',units='m',    &
                               FillValue=fv,missing_value=mv, &
                               netcdf_real=NF90_DOUBLE)
      end select
   end if

   if (init3d) then
      status = nf90_def_var(ncid,zname,NF90_DOUBLE, (/ z_dim /),id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,               &
                                    "init_grid_ncdf()","z -")
      call set_attributes(ncid,id,units=zunits)
   endif


!  bathymetry
   status = nf90_def_var(ncid,'bathymetry',NF90_DOUBLE,f2_dims,id)
   if (status .ne. NF90_NOERR) call netcdf_error(status,                  &
                                  "init_grid_ncdf()","bathymetry -")
   fv = h_missing
   mv = h_missing
   vr(1) = -5.
   vr(2) = 4000.
   call set_attributes(ncid,id,                                         &
                       long_name='bathymetry',units='m',                &
                       netcdf_real=NF90_DOUBLE, &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   if (save_masks) then

      status = nf90_def_var(ncid,'t_mask',NF90_INT,f2_dims,id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,             &
                                  "init_grid_ncdf()","t_mask")
      call set_attributes(ncid,id,long_name='mask for T-points')

      status = nf90_def_var(ncid,'u_mask',NF90_INT,f2_dims,id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,             &
                                  "init_grid_ncdf()","u_mask")
      call set_attributes(ncid,id,long_name='mask for U-points')

      status = nf90_def_var(ncid,'v_mask',NF90_INT,f2_dims,id)
      if (status .ne. NF90_NOERR) call netcdf_error(status,             &
                                  "init_grid_ncdf()","v_mask")
      call set_attributes(ncid,id,long_name='mask for V-points')

   end if


   return
   end subroutine init_grid_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
