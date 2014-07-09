#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise 3D netCDF variables
!
! !INTERFACE:
   subroutine init_3d_ncdf(fn,title,starttime,runtype)
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use exceptions
   use ncdf_common
   use ncdf_3d
   use domain, only: ioff,joff
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: vert_cord
   use m2d, only: no_2d
   use m3d, only: calc_temp,calc_salt
   use nonhydrostatic, only: nonhyd_iters,bnh_filter,bnh_weight,calc_hs2d,sbnh_filter
#ifdef SPM
   use suspended_matter, only: spm_save
#endif
#ifdef GETM_BIO
   use bio_var, only: numc,var_names,var_units,var_long
#endif
#ifdef _FABM_
   use getm_fabm, only: model,fabm_calc
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,title,starttime
   integer, intent(in)                 :: runtype
!
! !DEFINED PARAMETERS:
   logical,    parameter               :: init3d=.true.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: n,rc
   integer                   :: scalar(1),f3_dims(3),f4_dims(4)
   REALTYPE                  :: fv,mv,vr(2)
   character(len=80)         :: history,ts
!EOP
!-------------------------------------------------------------------------
!BOC
!  create netCDF file
   err = nf90_create(fn, NF90_CLOBBER, ncid)
   if (err .NE. NF90_NOERR) go to 10

!  initialize all time-independent, grid related variables
   call init_grid_ncdf(ncid,init3d,x_dim,y_dim,z_dim)

!  define unlimited dimension
   err = nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dim)
   if (err .NE. NF90_NOERR) go to 10

!  netCDF dimension vectors
   f3_dims(3)= time_dim
   f3_dims(2)= y_dim
   f3_dims(1)= x_dim

   f4_dims(4)= time_dim
   f4_dims(3)= z_dim
   f4_dims(2)= y_dim
   f4_dims(1)= x_dim

!  gobal settings
   history = 'GETM, ver. '//RELEASE
   ts = 'seconds since '//starttime

!  time
   err = nf90_def_var(ncid,'time',NF90_DOUBLE,time_dim,time_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,time_id,units=trim(ts),long_name='time')

!  elevation
   err = nf90_def_var(ncid,'elev',NCDF_FLOAT_PRECISION,f3_dims,elev_id)
   if (err .NE. NF90_NOERR) go to 10
   fv = elev_missing
   mv = elev_missing
   vr(1) = -15.
   vr(2) =  15.
   call set_attributes(ncid,elev_id,long_name='elevation',units='m', &
                       FillValue=fv,missing_value=mv,valid_range=vr)


   if (save_fluxes) then

      fv = vel_missing
      mv = vel_missing
      vr(1) = -3.
      vr(2) =  3.

      err = nf90_def_var(ncid,'fluxu_adv',NCDF_FLOAT_PRECISION,f3_dims,fluxu_adv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,fluxu_adv_id,long_name='avg. grid-related volume flux in local x-direction (U-point)',units='m3/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'fluxv_adv',NCDF_FLOAT_PRECISION,f3_dims,fluxv_adv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,fluxv_adv_id,long_name='avg. grid-related volume flux in local y-direction (V-point)',units='m3/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'fluxuu',NCDF_FLOAT_PRECISION,f4_dims,fluxuu_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,fluxuu_id,long_name='grid-related volume flux in local x-direction (U-point)',units='m3/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'fluxvv',NCDF_FLOAT_PRECISION,f4_dims,fluxvv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,fluxvv_id,long_name='grid-related volume flux in local y-direction (V-point)',units='m3/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'fluxw',NCDF_FLOAT_PRECISION,f4_dims,fluxw_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,fluxw_id,long_name='vertical volume flux (W-point)',units='m3/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

   end if


   fv = vel_missing
   mv = vel_missing
   vr(1) = -1.
   vr(2) =  1.

   if (save_vel2d) then

      err = nf90_def_var(ncid,'u_adv',NCDF_FLOAT_PRECISION,f3_dims,u_adv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,u_adv_id,long_name='avg. velocity in global x-direction (T-point)',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'v_adv',NCDF_FLOAT_PRECISION,f3_dims,v_adv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,v_adv_id,long_name='avg. velocity in global y-direction (T-point)',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

   end if

   if (save_vel3d) then

      err = nf90_def_var(ncid,'uu',NCDF_FLOAT_PRECISION,f4_dims,uu_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,uu_id,long_name='velocity in global x-direction (T-point)',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'vv',NCDF_FLOAT_PRECISION,f4_dims,vv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,vv_id,long_name='velocity in global y-direction (T-point)',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'w',NCDF_FLOAT_PRECISION,f4_dims,w_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,w_id,long_name='vertical velocity (T-point)',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

   if (save_taub) then

      !  bottom stress in x-direction
      err = nf90_def_var(ncid,'taubx',NCDF_FLOAT_PRECISION,f3_dims,taubx_id)
      if (err .NE. NF90_NOERR) go to 10
      fv = tau_missing
      mv = tau_missing
      vr(1) = -10.
      vr(2) =  10.
      call set_attributes(ncid,taubx_id,long_name='bottom stress (x)',units='Pa', &
           FillValue=fv,missing_value=mv,valid_range=vr)

      !  bottom stress in y-direction
      err = nf90_def_var(ncid,'tauby',NCDF_FLOAT_PRECISION,f3_dims,tauby_id)
      if (err .NE. NF90_NOERR) go to 10
      fv = tau_missing
      mv = tau_missing
      vr(1) = -10.
      vr(2) =  10.
      call set_attributes(ncid,tauby_id,long_name='bottom stress (y)',units='Pa', &
           FillValue=fv,missing_value=mv,valid_range=vr)
   endif

   if (save_h) then
      fv = hh_missing
      mv = hh_missing
      err = nf90_def_var(ncid,'h',NCDF_FLOAT_PRECISION,f4_dims,h_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,h_id,long_name='layer thickness',  &
                          units='m',FillValue=fv,missing_value=mv)
   end if


!  hydrostatic consistency criterion
   err = nf90_def_var(ncid,'hcc',NCDF_FLOAT_PRECISION,f4_dims(1:3),hcc_id)
   if (err .NE. NF90_NOERR) go to 10
   fv = -_ONE_
   mv = -_ONE_
   vr(1) = 0.
   vr(2) = 1.
   call set_attributes(ncid,hcc_id,  &
                       long_name='hcc',units=' ',          &
                       FillValue=fv,missing_value=mv,valid_range=vr)

#ifdef _MOMENTUM_TERMS_
   fv = vel_missing
   mv = vel_missing
   vr(1) = -3.
   vr(2) =  3.

   err = nf90_def_var(ncid,'tdv_u',NCDF_FLOAT_PRECISION,f4_dims,tdv_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,tdv_u_id,long_name='time tendency (u)',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'adv_u',NCDF_FLOAT_PRECISION,f4_dims,adv_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,adv_u_id,long_name='advection (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'vsd_u',NCDF_FLOAT_PRECISION,f4_dims,vsd_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,vsd_u_id,long_name='vertical stress divergence (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'hsd_u',NCDF_FLOAT_PRECISION,f4_dims,hsd_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,hsd_u_id,long_name='horizontal stress divergence (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'cor_u',NCDF_FLOAT_PRECISION,f4_dims,cor_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,cor_u_id,long_name='Coriolis term (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'epg_u',NCDF_FLOAT_PRECISION,f4_dims,epg_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,epg_u_id,long_name='extenal pressure gradient (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'ipg_u',NCDF_FLOAT_PRECISION,f4_dims,ipg_u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,ipg_u_id,long_name='internal pressure gradient (u).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'tdv_v',NCDF_FLOAT_PRECISION,f4_dims,tdv_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,tdv_v_id,long_name='time tendency (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'adv_v',NCDF_FLOAT_PRECISION,f4_dims,adv_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,adv_v_id,long_name='advection (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'vsd_v',NCDF_FLOAT_PRECISION,f4_dims,vsd_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,vsd_v_id,long_name='vertical stress divergence (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'hsd_v',NCDF_FLOAT_PRECISION,f4_dims,hsd_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,hsd_v_id,long_name='horizontal stress divergence (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'cor_v',NCDF_FLOAT_PRECISION,f4_dims,cor_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,cor_v_id,long_name='Coriolis term (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'epg_v',NCDF_FLOAT_PRECISION,f4_dims,epg_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,epg_v_id,long_name='external pressure gradient (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

   err = nf90_def_var(ncid,'ipg_v',NCDF_FLOAT_PRECISION,f4_dims,ipg_v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,ipg_v_id,long_name='internal pressure gradient (v).',units='m2/s2', &
                       FillValue=fv,missing_value=mv,valid_range=vr)
#endif

   if (save_s) then
      fv = salt_missing
      mv = salt_missing
      vr(1) =  0.
      vr(2) = 40.
      err = nf90_def_var(ncid,'salt',NCDF_FLOAT_PRECISION,f4_dims,salt_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,salt_id,long_name='salinity',units='PSU', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

   if (save_t) then
      fv = temp_missing
      mv = temp_missing
      vr(1) = -2.
      vr(2) = 40.
      err = nf90_def_var(ncid,'temp',NCDF_FLOAT_PRECISION,f4_dims,temp_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,temp_id,long_name='temperature',units='degC',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

   if (save_rho) then
      fv = rho_missing
      mv = rho_missing
      vr(1) =  0.
      vr(2) = 30.
      err = nf90_def_var(ncid,'sigma_t',NCDF_FLOAT_PRECISION,f4_dims,sigma_t_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,sigma_t_id,long_name='sigma_t',units='kg/m3',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

   if (save_strho) then
      if (save_rad) then
         fv = rad_missing
         mv = rad_missing
         vr(1) =  0.
         vr(2) = 1354.
         err = nf90_def_var(ncid,'radiation',NCDF_FLOAT_PRECISION,f4_dims,rad_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,rad_id,long_name='radiation',units='W/m2',&
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if
   end if

#ifndef NO_BAROCLINIC
   if (calc_stirr .and. save_stirr) then
      fv = stirr_missing
      mv = stirr_missing
      vr(1) = -500.
      vr(2) =  500.

      err = nf90_def_var(ncid,'diffxx',NCDF_FLOAT_PRECISION,f4_dims,diffxx_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,diffxx_id,long_name='zonal stirring diffusivity',units='m2/s',&
                          FillValue=fv,missing_value=mv,valid_range=vr)

#ifndef SLICE_MODEL
      err = nf90_def_var(ncid,'diffyy',NCDF_FLOAT_PRECISION,f4_dims,diffyy_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,diffyy_id,long_name='meridional stirring diffusivity',units='m2/s',&
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'diffxy',NCDF_FLOAT_PRECISION,f4_dims,diffxy_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,diffxy_id,long_name='cross stirring diffusivity',units='m2/s',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
#endif
   end if
#endif

   if (save_turb) then

      if (save_tke) then
         fv = tke_missing
         mv = tke_missing
         vr(1) = 0.
         vr(2) = 0.2
         err = nf90_def_var(ncid,'tke',NCDF_FLOAT_PRECISION,f4_dims,tke_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,tke_id,long_name='TKE',units='m2/s2', &
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (save_num) then
         fv = num_missing
         mv = num_missing
         vr(1) = 0.
         vr(2) = 0.2
         err = nf90_def_var(ncid,'num',NCDF_FLOAT_PRECISION,f4_dims,num_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,num_id,long_name='viscosity',units='m2/s', &
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (save_nuh) then
         fv = nuh_missing
         mv = nuh_missing
         vr(1) = 0.
         vr(2) = 0.2
         err = nf90_def_var(ncid,'nuh',NCDF_FLOAT_PRECISION,f4_dims,nuh_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nuh_id,long_name='diffusivity',units='m2/s', &
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (save_eps) then
         fv = eps_missing
         mv = eps_missing
         vr(1) = 0.
         vr(2) = 0.2
         err = nf90_def_var(ncid,'diss',NCDF_FLOAT_PRECISION,f4_dims,eps_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,eps_id,long_name='dissipation',units='m2/s3',&
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if
   end if

   if (save_SS_NN) then

      fv = SS_missing
      mv = SS_missing
      vr(1) = 0.
      vr(2) = 0.01
      err = nf90_def_var(ncid,'SS',NCDF_FLOAT_PRECISION,f4_dims,SS_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,SS_id,long_name='shear frequency squared',units='s-2',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
#ifndef NO_BAROCLINIC
      fv = NN_missing
      mv = NN_missing
      vr(1) = -0.001
      vr(2) = 0.01
      err = nf90_def_var(ncid,'NN',NCDF_FLOAT_PRECISION,f4_dims,NN_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,NN_id,long_name='buoyancy frequency squared', &
                          units='s-2',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
#endif

   end if

   if (do_numerical_analyses_3d) then

      fv = nummix_missing
      mv = nummix_missing
      vr(1) = -100.0
      vr(2) = 100.0
      err = nf90_def_var(ncid,'numdis_3d',NCDF_FLOAT_PRECISION,f4_dims,nd3d_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,nd3d_id, &
          long_name='numerical dissipation', &
          units='W/kg',&
          FillValue=fv,missing_value=mv,valid_range=vr)

#ifdef _NUMERICAL_ANALYSES_OLD_
      err = nf90_def_var(ncid,'numdis_3d_old',NCDF_FLOAT_PRECISION,f4_dims,nd3do_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,nd3do_id, &
          long_name='numerical dissipation (old)', &
          units='W/kg',&
          FillValue=fv,missing_value=mv,valid_range=vr)
#endif

      err = nf90_def_var(ncid,'phydis_3d',NCDF_FLOAT_PRECISION,f4_dims,pd3d_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,pd3d_id, &
          long_name='physical dissipation', &
          units='W/kg',&
          FillValue=fv,missing_value=mv,valid_range=vr)

      if (calc_salt) then
         err = nf90_def_var(ncid,'nummix_S',NCDF_FLOAT_PRECISION,f4_dims,nmS_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nmS_id, &
             long_name='numerical mixing of salinity', &
             units='psu**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)

#ifdef _NUMERICAL_ANALYSES_OLD_
         err = nf90_def_var(ncid,'nummix_S_old',NCDF_FLOAT_PRECISION,f4_dims,nmSo_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nmSo_id, &
             long_name='numerical mixing of salinity (old)', &
             units='psu**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
#endif

         err = nf90_def_var(ncid,'phymix_S',NCDF_FLOAT_PRECISION,f4_dims,pmS_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,pmS_id, &
             long_name='physical mixing of salinity', &
             units='psu**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (calc_temp) then
         err = nf90_def_var(ncid,'nummix_T',NCDF_FLOAT_PRECISION,f4_dims,nmT_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nmT_id, &
             long_name='numerical mixing of temperature', &
             units='degC**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)

#ifdef _NUMERICAL_ANALYSES_OLD_
         err = nf90_def_var(ncid,'nummix_T_old',NCDF_FLOAT_PRECISION,f4_dims,nmTo_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nmTo_id, &
             long_name='numerical mixing of temperature (old)', &
             units='degC**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
#endif

         err = nf90_def_var(ncid,'phymix_T',NCDF_FLOAT_PRECISION,f4_dims,pmT_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,pmT_id, &
             long_name='physical mixing of temperature', &
             units='degC**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

   end if

   if (nonhyd_method .ne. 0) then
      fv = bnh_missing
      mv = bnh_missing
      if (runtype.eq.2 .or. nonhyd_method.eq.1) then
         vr(1) = -10.
         vr(2) = 10.
         err = nf90_def_var(ncid,'bnh',NCDF_FLOAT_PRECISION,f4_dims,bnh_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,bnh_id,long_name='nh buoyancy correction',units='m/s2',&
                             FillValue=fv,missing_value=mv,valid_range=vr)
         if (nonhyd_method .eq. 1) then
            err = nf90_put_att(ncid,bnh_id,'nonhyd_iters',nonhyd_iters)
            err = nf90_put_att(ncid,bnh_id,'bnh_filter',bnh_filter)
            if (bnh_filter .eq. 1 .or. bnh_filter .eq. 3) then
               err = nf90_put_att(ncid,bnh_id,'bnh_weight',bnh_weight)
            end if
            if (.not. no_2d) then
               if (calc_hs2d) then
                  err = nf90_put_att(ncid,bnh_id,'calc_hs2d','true')
               else
                  err = nf90_put_att(ncid,bnh_id,'calc_hs2d','false')
                  if (sbnh_filter) then
                     err = nf90_put_att(ncid,bnh_id,'sbnh_filter','true')
                  else
                     err = nf90_put_att(ncid,bnh_id,'sbnh_filter','false')
                  end if
               end if
            end if
         end if
      else
         vr(1) = -10./SMALL
         vr(2) =  10./SMALL
         err = nf90_def_var(ncid,'nhsp',NCDF_FLOAT_PRECISION,f4_dims,bnh_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,bnh_id,long_name='nh screening parameter',units=' ',&
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if
   end if

   if (Am_method.eq.AM_LES .and. save_Am_3d) then
      fv = Am_3d_missing
      mv = Am_3d_missing
      vr(1) = 0.
      vr(2) = 500.
      err = nf90_def_var(ncid,'Am_3d',NCDF_FLOAT_PRECISION,f4_dims,Am_3d_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,Am_3d_id,long_name='hor. eddy viscosity',units='m2/s',&
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

#ifdef SPM
   if (spm_save) then
      fv = spm_missing
      mv = spm_missing
      err = nf90_def_var(ncid,'spm_pool',NCDF_FLOAT_PRECISION,f3_dims,spmpool_id)
      if (err .NE. NF90_NOERR) go to 10
      vr(1) = 0.
      vr(2) = 10.
      call set_attributes(ncid,spmpool_id,long_name='bottom spm pool', &
                          units='kg/m2', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
      vr(1) =  0.
      vr(2) = 30.
      err = nf90_def_var(ncid,'spm',NCDF_FLOAT_PRECISION,f4_dims,spm_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,spm_id,  &
                          long_name='suspended particulate matter', &
                          units='kg/m3', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if
#endif

#ifdef GETM_BIO
   allocate(bio_ids(numc),stat=rc)
   if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (bio_ids)'
   STDERR numc
   fv = bio_missing
   mv = bio_missing
   vr(1) = _ZERO_
   vr(2) = 9999.
   do n=1,numc
      err = nf90_def_var(ncid,var_names(n),NCDF_FLOAT_PRECISION,f4_dims,bio_ids(n))
      if (err .NE.  NF90_NOERR) go to 10
      call set_attributes(ncid,bio_ids(n), &
                          long_name=trim(var_long(n)), &
                          units=trim(var_units(n)), &
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end do
#endif

#ifdef _FABM_
   if (fabm_calc) then
      allocate(fabm_ids(size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (fabm_ids)'
      do n=1,size(model%state_variables)
         err = nf90_def_var(ncid,model%state_variables(n)%name,NCDF_FLOAT_PRECISION,f4_dims,fabm_ids(n))
         if (err .NE.  NF90_NOERR) go to 10
         call set_attributes(ncid,fabm_ids(n), &
                          long_name    =trim(model%state_variables(n)%long_name), &
                          units        =trim(model%state_variables(n)%units),    &
                          FillValue    =model%state_variables(n)%missing_value,  &
                          missing_value=model%state_variables(n)%missing_value,  &
                          valid_min    =model%state_variables(n)%minimum,        &
                          valid_max    =model%state_variables(n)%maximum)
      end do

      allocate(fabm_ids_ben(size(model%bottom_state_variables)),stat=rc)
      if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (fabm_ids_ben)'
      do n=1,size(model%bottom_state_variables)
         err = nf90_def_var(ncid,model%bottom_state_variables(n)%name,NCDF_FLOAT_PRECISION,f3_dims,fabm_ids_ben(n))
         if (err .NE.  NF90_NOERR) go to 10
         call set_attributes(ncid,fabm_ids_ben(n), &
                          long_name    =trim(model%bottom_state_variables(n)%long_name), &
                          units        =trim(model%bottom_state_variables(n)%units),    &
                          FillValue    =model%bottom_state_variables(n)%missing_value,  &
                          missing_value=model%bottom_state_variables(n)%missing_value,  &
                          valid_min    =model%bottom_state_variables(n)%minimum,        &
                          valid_max    =model%bottom_state_variables(n)%maximum)
      end do

      allocate(fabm_ids_diag(size(model%diagnostic_variables)),stat=rc)
      if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (fabm_ids_diag)'
      do n=1,size(model%diagnostic_variables)
         err = nf90_def_var(ncid,model%diagnostic_variables(n)%name,NCDF_FLOAT_PRECISION,f4_dims,fabm_ids_diag(n))
         if (err .NE.  NF90_NOERR) go to 10
         call set_attributes(ncid,fabm_ids_diag(n), &
                          long_name    =trim(model%diagnostic_variables(n)%long_name), &
                          units        =trim(model%diagnostic_variables(n)%units),    &
                          FillValue    =model%diagnostic_variables(n)%missing_value,  &
                          missing_value=model%diagnostic_variables(n)%missing_value,  &
                          valid_min    =model%diagnostic_variables(n)%minimum,        &
                          valid_max    =model%diagnostic_variables(n)%maximum)
      end do

      allocate(fabm_ids_diag_hz(size(model%horizontal_diagnostic_variables)),stat=rc)
      if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (fabm_ids_diag_hz)'
      do n=1,size(model%horizontal_diagnostic_variables)
         err = nf90_def_var(ncid,model%horizontal_diagnostic_variables(n)%name,NCDF_FLOAT_PRECISION,f3_dims,fabm_ids_diag_hz(n))
         if (err .NE.  NF90_NOERR) go to 10
         call set_attributes(ncid,fabm_ids_diag_hz(n), &
                          long_name    =trim(model%horizontal_diagnostic_variables(n)%long_name), &
                          units        =trim(model%horizontal_diagnostic_variables(n)%units),    &
                          FillValue    =model%horizontal_diagnostic_variables(n)%missing_value,  &
                          missing_value=model%horizontal_diagnostic_variables(n)%missing_value,  &
                          valid_min    =model%horizontal_diagnostic_variables(n)%minimum,        &
                          valid_max    =model%horizontal_diagnostic_variables(n)%maximum)
      end do

      if (do_numerical_analyses_3d) then

         allocate(pmpel_ids(size(model%state_variables)),stat=rc)
         if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (pmpel_ids)'
         do n=1,size(model%state_variables)
            err = nf90_def_var(ncid,'phymix_'//trim(model%state_variables(n)%name),NCDF_FLOAT_PRECISION,f4_dims,pmpel_ids(n))
            if (err .NE.  NF90_NOERR) go to 10
            call set_attributes(ncid,pmpel_ids(n),                                                                &
                             long_name    ='physical mixing of '//trim(model%state_variables(n)%long_name), &
                             units        =trim(model%state_variables(n)%units//'**2/s'),                    &
                             FillValue    =model%state_variables(n)%missing_value,                           &
                             missing_value=model%state_variables(n)%missing_value)
         end do

         allocate(nmpel_ids(size(model%state_variables)),stat=rc)
         if (rc /= 0) stop 'init_3d_ncdf(): Error allocating memory (nmpel_ids)'
         do n=1,size(model%state_variables)
            err = nf90_def_var(ncid,'nummix_'//trim(model%state_variables(n)%name),NCDF_FLOAT_PRECISION,f4_dims,nmpel_ids(n))
            if (err .NE.  NF90_NOERR) go to 10
            call set_attributes(ncid,nmpel_ids(n),                                                                &
                             long_name    ='numerical mixing of '//trim(model%state_variables(n)%long_name), &
                             units        =trim(model%state_variables(n)%units//'**2/s'),                    &
                             FillValue    =model%state_variables(n)%missing_value,                           &
                             missing_value=model%state_variables(n)%missing_value)
         end do
      end if

   end if
#endif

!  globals
   err = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_put_att(ncid,NF90_GLOBAL,'history',trim(history))
   if (err .NE. NF90_NOERR) go to 10

   history = GIT_REVISION
   err = nf90_put_att(ncid,NF90_GLOBAL,'git',trim(history))
   if (err .NE. NF90_NOERR) go to 10

   history = FORTRAN_VERSION
   err = nf90_put_att(ncid,NF90_GLOBAL,'compiler',trim(history))
   if (err .NE. NF90_NOERR) go to 10

   ! leave define mode
   err = nf90_enddef(ncid)
   if (err .NE. NF90_NOERR) go to 10

   return

   10 FATAL 'init_3d_ncdf: ',nf90_strerror(err)
   stop 'init_3d_ncdf'
   end subroutine init_3d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
