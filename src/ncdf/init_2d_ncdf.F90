#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise 2D netCDf variables
!
! !INTERFACE:
   subroutine init_2d_ncdf(fn,title,starttime)
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use exceptions
   use ncdf_common
   use ncdf_2d
   use domain, only: imin,imax,jmin,jmax
   use domain, only: ioff,joff
   use meteo,  only: metforcing,calc_met
   use meteo,  only: fwf_method
   use m2d,    only: residual

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,title,starttime

! !DEFINED PARAMETERS:
   logical,    parameter               :: init3d=.false.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: err
   integer                   :: scalar(1),f2_dims(2),f3_dims(3)
   REALTYPE                  :: fv,mv,vr(2)
   character(len=80)         :: history,ts
!EOP
!-------------------------------------------------------------------------
!BOC
!  create netCDF file
   err = nf90_create(fn, NF90_CLOBBER, ncid)
   if (err .NE. NF90_NOERR) go to 10

!  initialize all time-independent, grid related variables
   call init_grid_ncdf(ncid,init3d,x_dim,y_dim)

!  allocate workspace
   allocate(ws(E2DFIELD),stat=err)
   if (err .ne. 0) call getm_error("init_2d_ncdf()",               &
                                   "error allocating ws")

!  define unlimited dimension
   err = nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dim)
   if (err .NE. NF90_NOERR) go to 10

!  netCDF dimension vectors
   f3_dims(3)= time_dim
   f3_dims(2)= y_dim
   f3_dims(1)= x_dim

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
   call set_attributes(ncid,elev_id,long_name='elevation',units='meters', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

!  velocities
   fv = vel_missing
   mv = vel_missing
   vr(1) = -3.
   vr(2) =  3.
!  zonal velocity
   err = nf90_def_var(ncid,'u',NCDF_FLOAT_PRECISION,f3_dims,u_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,u_id,long_name='zonal velocity',units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

!  meridional velocity
   err = nf90_def_var(ncid,'v',NCDF_FLOAT_PRECISION,f3_dims,v_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,v_id,long_name='meridional velocity',units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)
#if defined(CURVILINEAR)
!  rotated zonal velocity
   err = nf90_def_var(ncid,'urot',NCDF_FLOAT_PRECISION,f3_dims,urot_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,urot_id,long_name='rot. zonal velocity', &
                       units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)

!  rotated meridional velocity
   err = nf90_def_var(ncid,'vrot',NCDF_FLOAT_PRECISION,f3_dims,vrot_id)
   if (err .NE. NF90_NOERR) go to 10
   call set_attributes(ncid,vrot_id,long_name='rot. meridional velocity', &
                       units='m/s', &
                       FillValue=fv,missing_value=mv,valid_range=vr)
#endif

!  meteorology
   if (metforcing .and. save_meteo) then
      if (calc_met) then
         fv = vel_missing; mv = vel_missing; vr(1) = -50.; vr(2) =  50.
         err = nf90_def_var(ncid,'u10',NCDF_FLOAT_PRECISION,f3_dims,u10_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,u10_id,long_name='U10',units='m/s', &
                             FillValue=fv,missing_value=mv,valid_range=vr)
         err = nf90_def_var(ncid,'v10',NCDF_FLOAT_PRECISION,f3_dims,v10_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,v10_id,long_name='V10',units='m/s', &
                             FillValue=fv,missing_value=mv,valid_range=vr)

         fv = airp_missing; mv = airp_missing;
         vr(1) = 90.e3; vr(2) = 110.e3
         err = nf90_def_var(ncid,'airp',NCDF_FLOAT_PRECISION,f3_dims,airp_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,airp_id,  &
                             long_name='air pressure',units='Pascal', &
                             FillValue=fv,missing_value=mv,valid_range=vr)

         fv = t2_missing; mv = t2_missing; vr(1) = 0; vr(2) = 325.
         err = nf90_def_var(ncid,'t2',NCDF_FLOAT_PRECISION,f3_dims,t2_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,t2_id,  &
                             long_name='temperature (2m)',units='Kelvin', &
                             FillValue=fv,missing_value=mv,valid_range=vr)

         fv = hum_missing; mv = hum_missing; vr(1) = 0; vr(2) = 100.
         err = nf90_def_var(ncid,'hum',NCDF_FLOAT_PRECISION,f3_dims,hum_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,hum_id,  &
                             long_name='humidity',units='kg/kg', &
                             FillValue=fv,missing_value=mv,valid_range=vr)

         fv = tcc_missing; mv = tcc_missing; vr(1) = 0.; vr(2) = 1.
         err = nf90_def_var(ncid,'tcc',NCDF_FLOAT_PRECISION,f3_dims,tcc_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,tcc_id,  &
                             long_name='total cloud cover',units='fraction', &
                             FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      fv = stress_missing; mv = stress_missing; vr(1) = -1; vr(2) = 1.
      err = nf90_def_var(ncid,'tausx',NCDF_FLOAT_PRECISION,f3_dims,tausx_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,tausx_id,  &
                          long_name='surface stress - x',units='N/m2', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'tausy',NCDF_FLOAT_PRECISION,f3_dims,tausy_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,tausy_id,  &
                          long_name='surface stress - y',units='N/m2', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      fv = swr_missing; mv = swr_missing; vr(1) = 0; vr(2) = 1500.
      err = nf90_def_var(ncid,'swr',NCDF_FLOAT_PRECISION,f3_dims,swr_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,swr_id,  &
                          long_name='short wave radiation',units='W/m2', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      fv = shf_missing; mv = shf_missing; vr(1) = -1000; vr(2) = 1000.
      err = nf90_def_var(ncid,'shf',NCDF_FLOAT_PRECISION,f3_dims,shf_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,shf_id,  &
                          long_name='surface heat fluxes',units='W/m2', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      if (fwf_method .ge. 2) then
         fv = evap_missing; mv = evap_missing; vr(1) = -1.0; vr(2) = 1.
         err = nf90_def_var(ncid,'evap',NCDF_FLOAT_PRECISION,f3_dims,evap_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,evap_id,  &
                            long_name='evaporation',units='m/s', &
                            FillValue=fv,missing_value=mv,valid_range=vr)
      end if

      if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
         fv = precip_missing; mv = precip_missing; vr(1) = -1.; vr(2) = 1.
         err = nf90_def_var(ncid,'precip',NCDF_FLOAT_PRECISION,f3_dims,precip_id)
         if (err .NE. NF90_NOERR) go to 10
         call set_attributes(ncid,precip_id,  &
                            long_name='precipitation',units='m/s', &
                            FillValue=fv,missing_value=mv,valid_range=vr)
      end if

   end if

   if (residual .gt. 0) then
!     residual currents - u and v
      fv = vel_missing; mv = vel_missing; vr(1) = -3.; vr(2) =  3.
      err = nf90_def_var(ncid,'res_u',NCDF_FLOAT_PRECISION,f3_dims,res_u_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,res_u_id,long_name='res. u',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)

      err = nf90_def_var(ncid,'res_v',NCDF_FLOAT_PRECISION,f3_dims,res_v_id)
      if (err .NE. NF90_NOERR) go to 10
      call set_attributes(ncid,res_v_id,long_name='res. v',units='m/s', &
                          FillValue=fv,missing_value=mv,valid_range=vr)
   end if

#ifdef USE_BREAKS
      err = nf90_def_var(ncid,'break_stat',NF90_INT,f3_dims(1:2),break_stat_id)
      if (err .ne. NF90_NOERR) call netcdf_error(err,                  &
                                  "init_2d_ncdf()","break_stat")
      call set_attributes(ncid,break_stat_id, &
                          long_name='stats (emergency breaks)')
#endif

!  globals
   err = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_put_att(ncid,NF90_GLOBAL,'version',trim(history))
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

   10 FATAL 'init_2d_ncdf: ',nf90_strerror(err)
   stop 'init_2d_ncdf'
   end subroutine init_2d_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
