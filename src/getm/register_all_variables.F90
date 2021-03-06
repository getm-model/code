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
   call register_2d_variables()
#ifndef NO_3D
   call register_3d_variables(runtype)
#endif
#ifdef _FABM_
   call register_fabm_variables(fm)
#endif
#if 0
   call register_diagnostic_variables()
#endif
   call register_processed_variables(fm,runtype)

   return
   end subroutine do_register_all_variables
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
! !ROUTINE: register_2d_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_2d_variables()
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   logical :: used
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_2d_variables()'
!D(E2DFIELD)
!DU,DV
!z(E2DFIELD)
!zo(E2DFIELD)
!U(E2DFIELD)
!V(E2DFIELD)
!UEx(E2DFIELD)
!VEx(E2DFIELD)
!fU(E2DFIELD)
!fV(E2DFIELD)
!ru(E2DFIELD)
!rv(E2DFIELD)
!Uint(E2DFIELD)
!Vint(E2DFIELD)
!Uinto(E2DFIELD)
!Vinto(E2DFIELD)
!res_du(E2DFIELD)
!res_u(E2DFIELD)
!res_dv(E2DFIELD)
!res_v(E2DFIELD)
!kbk
!SlUx(E2DFIELD)
!SlVx(E2DFIELD)
!Slru(E2DFIELD)
!Slrv(E2DFIELD)
!zub(E2DFIELD)
!zvb(E2DFIELD)
!zub0(E2DFIELD)
!zvb0(E2DFIELD)
!An(E2DFIELD)
!AnX(E2DFIELD)
!fwf(E2DFIELD)
!fwf_int(E2DFIELD)
!EWbdy(jmax),ENbdy(imax),EEbdy(jmax),ESbdy(imax)


!  category - 2d
   call fm%register('z', 'm', 'sea surface elevation', standard_name='sea surface elevation', fill_value=-9999.0_rk, data2d=z(_2D_W_), category="2d", part_of_state=.true.)
   call fm%register('zo', 'm', 'sea surface elevation', standard_name='sea surface elevation', fill_value=-9999.0_rk, data2d=zo(_2D_W_), category="2d", output_level=output_level_debug)
   call fm%register('D', 'm', 'water depth', standard_name='water depth', fill_value=-10009.0_rk, data2d=D(_2D_W_), category="2d", part_of_state=.true.)
   call fm%register('U', 'm2/s', 'transport in local x-direction', standard_name='', data2d=U(_2D_W_), category='2d', output_level=output_level_debug, part_of_state=.true.)
   call fm%register('V', 'm2/s', 'transport in local y-direction', standard_name='', data2d=V(_2D_W_), category='2d', output_level=output_level_debug, part_of_state=.true.)

   return
   end subroutine register_2d_variables
!EOC

#ifndef NO_3D
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_3d_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_3d_variables(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use variables_3d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_3d_variables()'

!:: kmin(I2DFIELD)
!:: kumin(I2DFIELD)
!:: kvmin(I2DFIELD)
!:: kmin_pmz(I2DFIELD)
!:: kumin_pmz(I2DFIELD)
!:: kvmin_pmz(I2DFIELD)

!:: uu(I3DFIELD)
!:: vv(I3DFIELD)
!:: ww(I3DFIELD)

#ifdef _MOMENTUM_TERMS_
!:: tdv_u(I3DFIELD)
!:: adv_u(I3DFIELD)
!:: vsd_u(I3DFIELD)
!:: hsd_u(I3DFIELD)
!:: cor_u(I3DFIELD)
!:: epg_u(I3DFIELD)
!:: ipg_u(I3DFIELD)

!:: tdv_v(I3DFIELD)
!:: adv_v(I3DFIELD)
!:: vsd_v(I3DFIELD)
!:: hsd_v(I3DFIELD)
!:: cor_v(I3DFIELD)
!:: epg_v(I3DFIELD)
!:: ipg_v(I3DFIELD)
#endif
#ifdef STRUCTURE_FRICTION
!:: sf(I3DFIELD)
#endif
!:: ho(I3DFIELD)
!:: hn(I3DFIELD)
!:: huo(I3DFIELD)
!:: hun(I3DFIELD)
!:: hvo(I3DFIELD)
!:: hvn(I3DFIELD)
!:: hcc(I3DFIELD)
!:: uuEx(I3DFIELD)
!:: vvEx(I3DFIELD)
!:: num(I3DFIELD)
!:: nuh(I3DFIELD)

! 3D turbulent fields
!:: tke(I3DFIELD)
!:: eps(I3DFIELD)
!:: SS(I3DFIELD)
#ifndef NO_BAROCLINIC
! 3D baroclinic fields
!:: NN(I3DFIELD)
!:: S(I3DFIELD)
!:: T(I3DFIELD)
!:: rho(I3DFIELD)
!:: rad(I3DFIELD)
!:: buoy(I3DFIELD)
!:: alpha(I3DFIELD)
!:: beta(I3DFIELD)
!:: idpdx(I3DFIELD)
!:: idpdy(I3DFIELD)
!:: light(I3DFIELD)
#endif

#ifdef SPM
! suspended matter
!:: spm(I3DFIELD)
!:: spm_ws(I3DFIELD)
!:: spm_pool(I2DFIELD)
#endif

! 2D fields in 3D domain
!:: sseo(I2DFIELD)
!:: ssen(I2DFIELD)
!:: Dn(I2DFIELD)
!:: ssuo(I2DFIELD)
!:: ssun(I2DFIELD)
!:: ssvo(I2DFIELD)
!:: ssvn(I2DFIELD)
!:: Dun,Dvn

! 3D friction in 3D domain
!:: rru(I2DFIELD)
!:: rrv(I2DFIELD)
!:: taus(I2DFIELD)
!:: taubx(I2DFIELD)
!:: tauby(I2DFIELD)
!:: taub(I2DFIELD)

! light attenuation
!:: A(I2DFIELD)
!:: g1(I2DFIELD)
!:: g2(I2DFIELD)

!  category - 3d
   if (runtype .ge. 2) then
      call fm%register('zc', 'm', 'center coordinate', standard_name='', dimensions=(/id_dim_z/),data3d=zc(_3D_W_), category='grid', part_of_state=.false.)
      call fm%register('hn', 'm', 'layer thickness', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hn(_3D_W_), category='grid', part_of_state=.true.)
      call fm%register('hun', 'm', 'layer thickness - U-points', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hun(_3D_W_), category='grid', output_level=output_level_debug)
      call fm%register('hvn', 'm', 'layer thickness - V-points', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hvn(_3D_W_), category='grid', output_level=output_level_debug)
      call fm%register('ho', 'm', 'old layer thickness', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=ho(_3D_W_), category='grid', output_level=output_level_debug)
      call fm%register('ssen', 'm', 'elevation at T-points (3D)', standard_name='', data2d=ssen(_2D_W_), category='3d', fill_value=-9999.0_rk, output_level=output_level_debug, part_of_state=.true.)
      call fm%register('ssun', 'm', 'elevation at U-points (3D)', standard_name='', data2d=ssun(_2D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
      call fm%register('ssvn', 'm', 'elevation at V-points (3D)', standard_name='', data2d=ssvn(_2D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
      call fm%register('sseo', 'm', 'old elevation at T-points (3D)', standard_name='', data2d=sseo(_2D_W_), category='3d', fill_value=-9999.0_rk, output_level=output_level_debug, part_of_state=.true.)
      call fm%register('uu', 'm2/s', 'transport in local x-direction (3D)', standard_name='', dimensions=(/id_dim_z/), data3d=uu(_3D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
      call fm%register('vv', 'm2/s', 'transport in local y-direction (3D)', standard_name='', dimensions=(/id_dim_z/), data3d=vv(_3D_W_), category='3d', output_level=output_level_debug, part_of_state=.true.)
      call fm%register('ww', 'm/s', 'grid-related vertical velocity', standard_name='', dimensions=(/id_dim_z/), data3d=ww(_3D_W_), category='3d', output_level=output_level_debug)
      call fm%register('SS', 's-2', 'shear frequency squared', standard_name='', dimensions=(/id_dim_z/), data3d=SS(_3D_W_), category='3d', output_level=output_level_debug)

   end if

!  category - turbulence
   if (runtype .ge. 2) then
      call fm%register('tke' , 'm2/s2', 'TKE'        , standard_name='', dimensions=(/id_dim_z/), data3d=tke(_3D_W_), category='turbulence', output_level=output_level_debug)
      call fm%register('diss', 'm2/s3', 'dissipation', standard_name='', dimensions=(/id_dim_z/), data3d=eps(_3D_W_), category='turbulence', output_level=output_level_debug)
      call fm%register('num' , 'm2/s' , 'viscosity'  , standard_name='', dimensions=(/id_dim_z/), data3d=num(_3D_W_), category='turbulence', output_level=output_level_debug)
      call fm%register('nuh' , 'm2/s' , 'diffusivity', standard_name='', dimensions=(/id_dim_z/), data3d=nuh(_3D_W_), category='turbulence', output_level=output_level_debug)
   end if

#ifndef NO_BAROCLINIC
!  category - baroclinic
   if (runtype .ge. 3) then
      call fm%register('temp', 'Celsius', 'temperature', standard_name='', dimensions=(/id_dim_z/), fill_value=-9999.0_rk, data3d=T  (_3D_W_), category='baroclinic', part_of_state=.true.)
      call fm%register('salt', '1e-3'   , 'salinity'   , standard_name='', dimensions=(/id_dim_z/), fill_value=-9999.0_rk, data3d=S  (_3D_W_), category='baroclinic', part_of_state=.true.)
      call fm%register('rho' , 'kg/m3'  , 'density'    , standard_name='', dimensions=(/id_dim_z/), fill_value=-9999.0_rk, data3d=rho(_3D_W_), category='baroclinic', output_level=output_level_debug)
      call fm%register('NN', 's-2', 'buoyancy frequency squared', standard_name='', dimensions=(/id_dim_z/), data3d=NN(_3D_W_), category='baroclinic', output_level=output_level_debug)
      call fm%register('idpdx', 'm2/s2', 'baroclinic pressure gradient - x', standard_name='', dimensions=(/id_dim_z/),data3d=idpdx(_3D_W_), category='baroclinic', output_level=output_level_debug)
#ifndef SLICE_MODEL
      call fm%register('idpdy', 'm2/s2', 'baroclinic pressure gradient - y', standard_name='', dimensions=(/id_dim_z/),data3d=idpdy(_3D_W_), category='baroclinic', output_level=output_level_debug)
#endif
   end if
#endif

   return
   end subroutine register_3d_variables
!EOC
#endif

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
