#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: register_all_variables
!
! !INTERFACE:
   module register_all_variables
!
! !DESCRIPTION:
!
! !USES:
   use field_manager
   use variables_2d, only: register_2d_variables
   use variables_3d, only: register_3d_variables
#ifdef _FABM_
   use getm_fabm, only: register_fabm_variables
#endif
   use output_processing, only: register_processed_variables, finalize_register_processed_variables
   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: init_register_all_variables
   public :: do_register_all_variables
   public :: finalize_register_all_variables
!
! !PUBLIC DATA MEMBERS:
   type (type_field_manager), public, target :: fm
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !PRIVATE DATA MEMBERS
   integer,parameter :: rk = kind(_ONE_)
   character(len=16)         :: xname=''
   character(len=16)         :: xlongname=''
   character(len=16)         :: xunits=''
   character(len=16)         :: yname=''
   character(len=16)         :: ylongname=''
   character(len=16)         :: yunits=''
   character(len=16)         :: zname='sigma'
   character(len=64)         :: zlongname='sigma'
   character(len=16)         :: zunits='sigma'
   character(len=16),parameter :: lonname     = 'lonc'
   character(len=16),parameter :: lonlongname = 'longitude'
   character(len=16),parameter :: lonunits    = 'degrees_east'
   character(len=16),parameter :: latname     = 'latc'
   character(len=16),parameter :: latlongname = 'latitude'
   character(len=16),parameter :: latunits    = 'degrees_north'
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_register_all_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine init_register_all_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use domain
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_register_all_variables()'

   select case (grid_type)
      case (1)
         xname     = 'xc'
         xlongname = 'x'
         xunits    = 'm'
         yname     = 'yc'
         ylongname = 'y'
         yunits    = 'm'
      case (2)
         xname     = lonname
         xlongname = lonlongname
         xunits    = lonunits
         yname     = latname
         ylongname = latlongname
         yunits    = latunits
      case (3,4)
         xname     = 'xic'
         xlongname = 'xic'
         yname     = 'etac'
         ylongname = 'etac'
   end select

#ifndef NO_3D
   if (runtype .ge. 2) then
      select case (vert_cord)
         case (1)
            zname     = 'sigma'
            zlongname = 'sigma layers'
            zunits    = 'sigma_level'
         case (2)
            zname     = 'z'
            zlongname = 'geopotential'
            zunits    = 'm'
         case (3,4,5)
            zname  = 'level'
            zlongname  = 'general vertical coordinates'
            zunits = 'level'
         case default
      end select
   end if
#endif

!  register - dimensions
   call fm%register_dimension(trim(xname),imax-imin+1,global_length=iextr,offset=ioff,id=id_dim_lon)
   call fm%register_dimension(trim(yname),jmax-jmin+1,global_length=jextr,offset=joff,id=id_dim_lat)
#ifndef NO_3D
   if (runtype .ge. 2) then
      call fm%register_dimension(trim(zname),kmax+1,global_length=kmax,offset=-1,id=id_dim_z)
   end if
#endif

   call fm%register_dimension('time',id=id_dim_time)
   call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat/),append_by_default=(/id_dim_time/))

   return
   end subroutine init_register_all_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_register_all_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine do_register_all_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d, only: register_2d_variables
#ifndef NO_3D
   use variables_3d, only: register_3d_variables
   use m3d, only: calc_temp,calc_salt
#endif
#ifdef _FABM_
   use getm_fabm, only: register_fabm_variables
#endif
   use output_processing, only: register_processed_variables
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'do_register_all_variables()'
   call register_domain_variables(runtype)
   call register_meteo_variables()
   call register_2d_variables(fm)
#ifndef NO_3D
   call register_3d_variables(fm,runtype,calc_temp,calc_salt)
#endif
#ifdef _FABM_
   call register_fabm_variables(fm)
#endif
#if 0
   call register_diagnostic_variables()
#endif
   call register_processed_variables(fm)

   return
   end subroutine do_register_all_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: finalize_register_all_variables() - send optional variables.
!
! !INTERFACE:
   subroutine finalize_register_all_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d, only: finalize_register_2d_variables
#ifndef NO_3D
   use variables_3d, only: finalize_register_3d_variables
   use m3d, only: calc_temp,calc_salt
#endif
#ifdef _FABM_
   use getm_fabm, only: finalize_register_fabm_variables
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'finalize_register_all_variables()'

   call finalize_register_2d_variables(fm)
#ifndef NO_3D
   call finalize_register_3d_variables(fm,calc_temp,calc_salt)
#endif
#ifdef _FABM_
   call finalize_register_fabm_variables(fm)
#endif

   return
   end subroutine finalize_register_all_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_domain_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_domain_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use domain
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_domain_variables()'

!  register - domain
   call fm%register(trim(xname),trim(xunits),trim(xlongname),dimensions=(/id_dim_lon/),no_default_dimensions=.true.,data1d=xcord(_IRANGE_NO_HALO_),coordinate_dimension=id_dim_lon,output_level=output_level_debug)
   call fm%register(trim(yname),trim(yunits),trim(ylongname),dimensions=(/id_dim_lat/),no_default_dimensions=.true.,data1d=ycord(_JRANGE_NO_HALO_),coordinate_dimension=id_dim_lat,output_level=output_level_debug)
#ifndef NO_3D
   if (runtype .ge. 2) then
      call fm%register(trim(zname),trim(zunits),trim(zlongname),dimensions=(/id_dim_z/),no_default_dimensions=.true.,data1d=ga,coordinate_dimension=id_dim_z,output_level=output_level_debug)
   end if
#endif

   if (have_lonlat .and. grid_type.ne.2) then
      call fm%register(trim(lonname),trim(lonunits),trim(lonlongname),dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=lonc(_2D_W_), category='domain',output_level=output_level_required)
      call fm%register(trim(latname),trim(latunits),trim(latlongname),dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=latc(_2D_W_), category='domain',output_level=output_level_required)
   end if

   call fm%register('bathymetry', 'm', 'bathymetry', standard_name='bathymetry', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., fill_value=-10._rk, data2d=H(_2D_W_), category='domain',output_level=output_level_required)

!  register -  metric
   call fm%register('dxc', 'm', 'dx at T-points', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., data2d=dxc(_2D_W_), category="metrics", output_level=output_level_debug)
   call fm%register('dyc', 'm', 'dy at T-points', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., data2d=dyc(_2D_W_), category="metrics", output_level=output_level_debug)

   return
   end subroutine register_domain_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_meteo_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_meteo_variables()
!
! !DESCRIPTION:
!
! !USES:
   use domain
   use meteo
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_meteo_variables()'

   if (metforcing) then
      if (calc_met) then
         call fm%register('airp', 'Pa', 'air pressure', standard_name='', data2d=airp(_2D_W_), category="meteo/in", output_level=output_level_debug)
         call fm%register('t2', 'Celcius', '2m air temperature', standard_name='', data2d=t2(_2D_W_), category="meteo/in", output_level=output_level_debug)
         call fm%register('u10', 'm/s', '10m wind (x)', standard_name='', data2d=u10(_2D_W_), category="meteo/in", output_level=output_level_debug)
         call fm%register('v10', 'm/s', '10m wind (y)', standard_name='', data2d=v10(_2D_W_), category="meteo/in", output_level=output_level_debug)
!:: hum
         call fm%register('tcc', '', 'total cloud cover', standard_name='', data2d=tcc(_2D_W_), category="meteo/in", output_level=output_level_debug)
         ! fwf_method = 2, 3 - precipitation read from file
         if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
            call fm%register('precip', 'm/s', 'precipitation', standard_name='', data2d=precip(_2D_W_), category="meteo/in", output_level=output_level_debug)
         end if
         ! fwf_method = 2 - evaporation read from file
         if (fwf_method .eq. 2) then
            call fm%register('evap', 'm/s', 'evaporation', standard_name='', data2d=evap(_2D_W_), category="meteo/in", output_level=output_level_debug)
         end if
      end if
      call fm%register('swr', 'W', 'short wave radiation', standard_name='', data2d=swr(_2D_W_), category="meteo/out", output_level=output_level_debug)
      call fm%register('shf', 'W', 'surface heat flux', standard_name='', data2d=shf(_2D_W_), category="meteo/out", output_level=output_level_debug)
      call fm%register('tausx', 'Pa', 'wind stress (x)', standard_name='', data2d=tausx(_2D_W_), category="meteo/out", output_level=output_level_debug)
      call fm%register('tausy', 'Pa', 'wind stress (y)', standard_name='', data2d=tausy(_2D_W_), category="meteo/out", output_level=output_level_debug)
      call fm%register('albedo', '', 'albedo', standard_name='', data2d=albedo(_2D_W_), category="meteo/out", output_level=output_level_debug)
      call fm%register('zenith_angle', 'degrees', 'solar zenith angle', standard_name='', data2d=zenith_angle(_2D_W_), category="meteo/out", output_level=output_level_debug)
      ! fwf_method = 3 - evaporation calculated
      if (fwf_method .eq. 3 .or. fwf_method .eq. 4) then
         call fm%register('evap', 'm/s', 'evaporation', standard_name='', data2d=evap(_2D_W_), category="meteo/out", output_level=output_level_debug)
      end if
   end if

!:: airp_old,airp_new
!:: tausx_old,tausy_old
!:: d_airp,d_tausx,d_tausy
!:: tcc_old,tcc_new
!:: swr_old,shf_old
!:: d_tcc,d_swr,d_shf
!:: evap_old,precip_old
!:: d_evap,d_precip

   return
   end subroutine register_meteo_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_diagnostic_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_diagnostic_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_diagnostic_variables() - non so-far'

   return
   end subroutine register_diagnostic_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: finalize_register_all_variables() - send optional variables.
!
! !INTERFACE:
   subroutine finalize_register_all_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d
#ifndef NO_3D
   use variables_3d
#endif
#ifdef _FABM_
   use getm_fabm
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'finalize_register_all_variables()'

   call finalize_register_processed_variables(fm)

   return
   end subroutine finalize_register_all_variables
!EOC

!-----------------------------------------------------------------------

   end module register_all_variables

!-----------------------------------------------------------------------
! Copyright (C) 2015 - Bolding & Bruggeman ApS                         !
!-----------------------------------------------------------------------
