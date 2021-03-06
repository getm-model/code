#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise 3D netCDF variables
!
! !INTERFACE:
   subroutine init_3d_ncdf(fn,title,starttime)
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
   use m3d, only: calc_temp,calc_salt
#ifdef SPM
   use suspended_matter, only: spm_save
#endif
#ifdef GETM_BIO
   use bio_var, only: numc,var_names,var_units,var_long
#endif
#ifdef _FABM_
   use getm_fabm, only: model,fabm_calc,output_none
#endif
   use getm_version
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,title,starttime
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
   history = 'GETM - www.getm.eu'
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

!  depth integrated zonal velocity
   err = nf90_def_var(ncid,'u',NCDF_FLOAT_PRECISION,f3_dims,u_id)
   if (err .NE. NF90_NOERR) go to 10
   fv = vel_missing
   mv = vel_missing
   vr(1) = -3.
   vr(2) =  3.
   call set_attributes(ncid,u_id,long_name='int. zonal vel.',units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

!  depth integrated meridional velocity
   err = nf90_def_var(ncid,'v',NCDF_FLOAT_PRECISION,f3_dims,v_id)
   if (err .NE. NF90_NOERR) go to 10
   fv = vel_missing
   mv = vel_missing
   vr(1) = -3.
   vr(2) =  3.
   call set_attributes(ncid,v_id,long_name='int. meridional vel.',units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)


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

   select case (vert_cord)
      case (_SIGMA_COORDS_)
      case (_Z_COORDS_)
         call getm_error("init_3d_ncdf()","saving of z-levels disabled")
      case (_GENERAL_COORDS_,_HYBRID_COORDS_,_ADAPTIVE_COORDS_)
         fv = hh_missing
         mv = hh_missing
         err = nf90_def_var(ncid,'h',NCDF_FLOAT_PRECISION,f4_dims,h_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,h_id,long_name='layer thickness',  &
                             units='m',FillValue=fv,missing_value=mv)
      case default
   end select


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

   if (save_vel) then
      fv = vel_missing
      mv = vel_missing
      vr(1) = -3.
      vr(2) =  3.

!     zonal velocity
      err = nf90_def_var(ncid,'uu',NCDF_FLOAT_PRECISION,f4_dims,uu_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,uu_id,long_name='zonal vel.',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

!     meridional velocity
      err = nf90_def_var(ncid,'vv',NCDF_FLOAT_PRECISION,f4_dims,vv_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,vv_id,long_name='meridional vel.',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

!     vertical velocity
      err = nf90_def_var(ncid,'w',NCDF_FLOAT_PRECISION,f4_dims,w_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,w_id,long_name='vertical vel.',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

#ifdef _MOMENTUM_TERMS_
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

#if defined(CURVILINEAR)
!     rotated zonal velocity
      err = nf90_def_var(ncid,'uurot',NCDF_FLOAT_PRECISION,f4_dims,uurot_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,uurot_id,long_name='rot. zonal vel.', &
                          units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

!     rotated meridional velocity
      err = nf90_def_var(ncid,'vvrot',NCDF_FLOAT_PRECISION,f4_dims,vvrot_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,vvrot_id,long_name='rot. meridional vel.', &
                          units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
#endif
   end if

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

   if (save_numerical_analyses) then

      fv = nummix_missing
      mv = nummix_missing
      vr(1) = -100.0
      vr(2) = 100.0
      err = nf90_def_var(ncid,'numdis3d',NCDF_FLOAT_PRECISION,f4_dims,nm3d_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,nm3d_id, &
          long_name='numerical dissipation', &
          units='W/kg',&
          FillValue=fv,missing_value=mv,valid_range=vr)

      if (calc_salt) then
         err = nf90_def_var(ncid,'nummix3d_S',NCDF_FLOAT_PRECISION,f4_dims,nm3dS_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nm3dS_id, &
             long_name='numerical mixing of salinity', &
             units='psu**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)

         err = nf90_def_var(ncid,'phymix3d_S',NCDF_FLOAT_PRECISION,f4_dims,pm3dS_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,pm3dS_id, &
             long_name='physical mixing of salinity', &
             units='psu**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (calc_temp) then
         err = nf90_def_var(ncid,'nummix3d_T',NCDF_FLOAT_PRECISION,f4_dims,nm3dT_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,nm3dT_id, &
             long_name='numerical mixing of temperature', &
             units='degC**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)

         err = nf90_def_var(ncid,'phymix3d_T',NCDF_FLOAT_PRECISION,f4_dims,pm3dT_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,pm3dT_id, &
             long_name='physical mixing of temperature', &
             units='degC**2/s',&
             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

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
      fabm_ids = -1
      do n=1,size(model%state_variables)
         if (model%state_variables(n)%output==output_none) cycle
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
      fabm_ids_ben = -1
      do n=1,size(model%bottom_state_variables)
         if (model%bottom_state_variables(n)%output==output_none) cycle
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
      fabm_ids_diag = -1
      do n=1,size(model%diagnostic_variables)
         if (model%diagnostic_variables(n)%output==output_none) cycle
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
      fabm_ids_diag_hz = -1
      do n=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(n)%output==output_none) cycle
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
   end if
#endif

!  globals
   err = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_put_att(ncid,NF90_GLOBAL,'model ',trim(history))
   if (err .NE. NF90_NOERR) go to 10

#if 0
   err = nf90_put_att(ncid,NF90_GLOBAL,'git hash:   ',trim(git_commit_id))
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_put_att(ncid,NF90_GLOBAL,'git branch: ',trim(git_branch_name))
   if (err .NE. NF90_NOERR) go to 10
#endif

!   history = FORTRAN_VERSION
!   err = nf90_put_att(ncid,NF90_GLOBAL,'compiler',trim(history))
!   if (err .NE. NF90_NOERR) go to 10

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
