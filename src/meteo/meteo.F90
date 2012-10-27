#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: meteo - provides meteorological forcing for \emph{getm}.
!
! !INTERFACE:
   module meteo
!
! !DESCRIPTION:
!  The meteo module provides meteorological forcing for \emph{getm}.
!  The main role of the module is to supply the following 4 fields with
!  sane values - \emph{airp, tausx, tausy, swr} and \emph{shf} - i.e.
!  air pressure [$Pa$], surface stresses [$N/m^2$], short wave radiation
!  [$W/m^2$] and surface heat fluxes [$W/m^2$] on the computational grid.
!  The module provides 3 public functions - \emph{init\_meteo()},
!  \emph{do\_meteo()} and \emph{clean\_meteo()} - and a number of public
!  data members to hold the actual meteorological fields. Also included
!  in the module are various constants related to calculating the
!  meteorological forcing.
!  Information about the calculation domain is obtained from the module
!  \emph{domain} and time related information comes from the module
!  \emph{time}.
!  The meteo module is initialised via a call to \emph{init\_meteo()}
!  that will read a namelist providing all necessary information. Memory
!  allocation is also done in \emph{init\_meteo()}.
!  Obtaining the actual forcing - either by reading from a file or
!  calculating is done via calls to \emph{do\_meteo()}. The actual
!  reading of external data from files is separated completely from
!  \emph{do\_meteo()} and is done in the main time loop via a call to
!  \emph{do\_input()} where all external file input is handled.
!  \emph{meteo} supplies 3 variables which can be used by routines for
!  reading variables. \emph{new\_meteo} is a logical switch which should
!  be set to .true. when new fields have been read. \emph{t\_1} and
!  \emph{t\_2} holds the time (in seconds) since the model run of the two
!  fields surrounding the actual model time - to be used by the temporal
!  interpolation. Finally \emph{clean\_meteo()} should be called when
!  the simulation is over as part of the overall procedure of finalising
!  the model run.
!
! !SEE ALSO:
!  short_wave_radiation.F90, fluxes.F90, exchange_coefficients.F90
!
! !USES:
   use time, only: yearday,secondsofday,timestep
   use halo_zones, only : H_TAG,update_2d_halo,wait_halo
   use domain, only: imin,imax,jmin,jmax,lonc,latc,convc,az
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_meteo, do_meteo, clean_meteo
!
! !PUBLIC DATA MEMBERS:
   character(LEN = PATH_MAX), public   :: meteo_file
   logical, public                     :: metforcing=.false.
   logical, public                     :: on_grid=.true.
   logical, public                     :: calc_met=.false.
   integer, public                     :: met_method
   integer, public                     :: fwf_method=0
   REALTYPE, public                    :: evap_factor = _ONE_
   REALTYPE, public                    :: precip_factor = _ONE_
   REALTYPE, public                    :: w,L,rho_air,qs,qa,ea,es
   REALTYPE, public, dimension(:,:), allocatable  :: u10,v10,t2,hum
   REALTYPE, public, dimension(:,:), pointer :: airp,tausx,tausy
   REALTYPE, public, dimension(:,:), pointer :: shf,swr=>null(),tcc
   REALTYPE, public, dimension(:,:), pointer :: evap,precip
   REALTYPE, public, dimension(:,:), pointer :: sst
   logical, public                           :: nudge_sst=.false.
   REALTYPE, public                          :: sst_const=-_ONE_
   REALTYPE, public                    :: cd_mom,cd_heat,cd_latent
   REALTYPE, public                    :: cd_precip = _ZERO_
   REALTYPE, public                    :: t_1=-_ONE_,t_2=-_ONE_
   logical, public                     :: new_meteo=.false.
   integer, public                     :: hum_method=-1
!
! !DEFINED PARAMETERS:
   REALTYPE,public,parameter           :: cpa=1008.d0  !AS that is not exact !
   REALTYPE,public,parameter           :: KELVIN=273.15d0
   REALTYPE,public,parameter           :: emiss=0.97d0
   REALTYPE,public,parameter           :: bolz=5.67d-8
!  REALTYPE,public,parameter           :: cpa=1004.67 ! specific heat of dry air- correct
   REALTYPE,public,parameter           :: cpw=4192.d0   ! specific heat of sea water
   REALTYPE,public,parameter           :: rho_precip = 1000.0d0
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: meteo_ramp=0,metfmt=2
   REALTYPE                  :: tx= _ZERO_ ,ty= _ZERO_
   REALTYPE                  :: swr_const= _ZERO_ ,shf_const= _ZERO_
   REALTYPE                  :: evap_const= _ZERO_ ,precip_const= _ZERO_
   REALTYPE, dimension(:,:), allocatable :: tausx_const,tausy_const
   REALTYPE, dimension(:,:), pointer     :: airp_new,d_airp
   REALTYPE, dimension(:,:), pointer     :: tausx_new,d_tausx
   REALTYPE, dimension(:,:), pointer     :: tausy_new,d_tausy
   REALTYPE, dimension(:,:), pointer     :: shf_new,d_shf
   REALTYPE, dimension(:,:), pointer     :: swr_new,d_swr
   REALTYPE, dimension(:,:), pointer     :: tcc_new,d_tcc
   REALTYPE, dimension(:,:), pointer     :: evap_new,d_evap
   REALTYPE, dimension(:,:), pointer     :: precip_new,d_precip
   REALTYPE, dimension(:,:), pointer     :: sst_new,d_sst
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_meteo - initialise the \emph{meteo} module.
!
! !INTERFACE:
   subroutine init_meteo(hotstart_method)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Basically reads the namelist \emph{meteo} from unit NAMLST.  According to
!  the content of the namelist various variables are allocated and initialised.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: hotstart_method
!
! !REVISION HISTORY:
!
!  See log for module.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,rc
   REALTYPE                  :: sinconv,cosconv
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029d0
   REALTYPE, parameter       :: deg2rad=pi/180
   namelist /meteo/ metforcing,on_grid,calc_met,met_method,fwf_method, &
                    meteo_ramp,metfmt,meteo_file, &
                    tx,ty,swr_const,shf_const,evap_const,precip_const, &
                    sst_const,precip_factor,evap_factor
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_meteo() # ',Ncall
#endif
   LEVEL1 'init_meteo'
   read(NAMLST,meteo)

   LEVEL2 'Metforcing=',metforcing
   if (metforcing) then
      select case (met_method)
            case (1)
               LEVEL2 'Constant forcing is used:'
               LEVEL3 'tx     = ',tx
               LEVEL3 'ty     = ',ty
               LEVEL3 'swr    = ',swr_const
               LEVEL3 'shf    = ',shf_const
            case (2)
               if(on_grid) then
                  LEVEL2 'Meteorological fields are on the computational grid'
               else
                  LEVEL2 'Meteorological fields needs to be interpolated'
               end if
               if(calc_met) then
                  LEVEL2 'Stresses and fluxes will be calculated'
               else
                  LEVEL2 'Stresses and fluxes are already calculated'
               end if
         case default
      end select

      select case (fwf_method)
            case (1)
               LEVEL2 'Constant evaporation/precipitation'
               LEVEL3 'evap   = ',evap_const
               LEVEL3 'precip = ',precip_const
            case (2)
               LEVEL2 'Evaporation/precipitation read from file'
            case (3)
               LEVEL2 'Precipitation read from file'
               LEVEL2 'Evaporation calculated'
            case (4)
               LEVEL2 'No precipitation'
               LEVEL2 'Evaporation calculated'
         case default
      end select

      if (meteo_ramp .gt. 1) then
         LEVEL2 'meteo_ramp=',meteo_ramp
         select case(hotstart_method)
            case (1)
               LEVEL3 'WARNING: re-start ramp for meteo'
            case (2)
               LEVEL3 'WARNING: no re-start of ramp for meteo'
         end select
      end if
   end if

! Allocates memory for the public data members

   allocate(airp(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (airp)'
   airp = _ZERO_

   allocate(tausx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (tausx)'
   allocate(tausy(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (tausy)'
   allocate(shf(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (shf)'
   allocate(swr(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (swr)'

   if (met_method .eq. 1) then
!     Rotation of wind stress due to grid convergence
      allocate(tausx_const(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausx_const)'
      allocate(tausy_const(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausy_const)'
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            sinconv=sin(-convc(i,j)*deg2rad)
            cosconv=cos(-convc(i,j)*deg2rad)
            tausx_const(i,j)= tx*cosconv+ty*sinconv
            tausy_const(i,j)=-tx*sinconv+ty*cosconv
         end do
      end do
      tausx = tausx_const
      tausy = tausy_const
      shf   = shf_const
      swr   = swr_const
   else
      tausx = _ZERO_
      tausy = _ZERO_
      shf   = _ZERO_
      swr   = _ZERO_
   end if

   allocate(evap(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (evap)'
   allocate(precip(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (precip)'
   if (fwf_method .eq. 1) then
      evap = evap_const
      precip = precip_const
   else
      evap = _ZERO_
      precip = _ZERO_
   end if


   if (metforcing .and. met_method.eq.2) then

      if (calc_met) then
         allocate(u10(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (u10)'
         u10 = _ZERO_

         allocate(v10(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (v10)'
         v10 = _ZERO_

         allocate(t2(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (t2)'
         t2 = _ZERO_

         allocate(hum(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (hum)'
         hum = _ZERO_

         allocate(tcc(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (tcc)'
         tcc = _ZERO_
      end if

      allocate(airp_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (airp_new)'
      allocate(d_airp(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_airp)'

      allocate(tausx_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausx_new)'
      allocate(d_tausx(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tausx)'

      allocate(tausy_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausy_new)'
      allocate(d_tausy(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tausy)'

      allocate(shf_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (shf_new)'
      allocate(d_shf(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_shf)'

      if (calc_met) then
         allocate(tcc_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (tcc_new)'
         allocate(d_tcc(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tcc)'
      else
         allocate(swr_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (swr_new)'
         allocate(d_swr(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_swr)'
      end if

      if (fwf_method .ge. 2) then
         allocate(evap_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (evap_new)'
         allocate(d_evap(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_evap)'
      end if

      if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
         allocate(precip_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (precip_new)'
         allocate(d_precip(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_precip)'
      end if

   end if
#ifdef DEBUG
   write(debug,*) 'Leaving init_meteo()'
   write(debug,*)
#endif
   return
   end subroutine init_meteo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_meteo - update the meteo forcing
!
! !INTERFACE:
   subroutine do_meteo(n,sst_model)
!$ use omp_lib
!
! !DESCRIPTION:
!  Should be called once every time step to update the meteorological forcing.
!  \emph{do\_meteo()} is called with two arguments - n is the loop number and
!  the sea surface temperature.
!  The modelled SST \emph{sst_model} is only used in the case where fluxes and stresses are
!  calculated as part of the model simulation.
!  The forcing can be obtained in 3 different way - using constant values,
!  using pre-calculated stresses and heat-fluxes or by calculating the
!  fluxes as part of the model integration. In all 3 cases the following
!  fields are the result \emph{, airp, tausx, tausy, swr} and \emph{shf} - i.e.
!  air pressure, stresses in x and y direction, short wave radiation and
!  surface heat fluxes.
!  The surface heat flux is the sum of the latent and sensible heat + the
!  net back radiation. Additional if available the surface freshwater fluxes
!  can be set const or read in from the meteo file. The unit in the meteo file
!  is assumed to be meter (/day).
!  The structure of this routine looks at first glance a bit more complicated
!  than should be necessary. The main reason is we need two fields in order to
!  do any time interpolation - which explains the use of \emph{first}.
!  In addition checks of the logical \emph{new\_meteo} is checked - set by the
!  reading subroutine.
!  Temporal interpolation is done in the principal variables.
!  It is possible to specify a soft start - via the meteo_ramp variable -
!  which is used to calculate a ramp (linearly from 0 to one over meteo_ramp
!  time steps).
!  To implement an use a different set of formulae for flux calculations
!  should be a matter of only changing the involved subroutines.
!
   use getm_timers, only: tic, toc, TIM_METEO
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in)                 :: n
   REALTYPE, optional, intent(inout)   :: sst_model(I2DFIELD)
!
! !REVISION HISTORY:
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,rc
   REALTYPE                  :: ramp,hh,t,t_minus_t2
   REALTYPE, save            :: deltm1
   REALTYPE                  :: short_wave_radiation
   logical,save              :: first=.true.
   REALTYPE, dimension(:,:), pointer :: airp_old,tausx_old,tausy_old
   REALTYPE, dimension(:,:), pointer :: shf_old,swr_old,tcc_old
   REALTYPE, dimension(:,:), pointer :: evap_old,precip_old
   REALTYPE, dimension(:,:), pointer :: sst_old
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_meteo() # ',Ncall
#endif
   call tic(TIM_METEO)
! OMP-NOTE: In this routine some loops, which have to do with read-in of
!    new meteo have not been threaded due to use of heap-allocated scalars
!    in exchange_coefficients() and fluxes().
!    However, for the cases I have tested, calls to short_wave_radiation
!    is by far the most expensive of the present routine.
!       BJB 2009-09-30.
   if (metforcing) then

      if (first) then
         if (nudge_sst) then
            if (met_method .eq. 2) then
               allocate(sst_new(E2DFIELD),stat=rc)
               if (rc /= 0) stop 'do_meteo: Error allocating memory (sst_new)'
               allocate(d_sst(E2DFIELD),stat=rc)
               if (rc /= 0) stop 'do_meteo: Error allocating memory (d_sst)'
            end if
         end if
      end if

      t = n*timestep
      hh = secondsofday*(_ONE_/3600)
      if (new_meteo) then
         if (.not. first) then
            deltm1 = _ONE_ / (t_2 - t_1)
         end if
      end if
      if (.not. first) then
         t_minus_t2 = t - t_2
      end if

      if(meteo_ramp.gt.1 .and. n.lt.meteo_ramp) then
         ramp = _ONE_*n/meteo_ramp
      else
         ramp = _ONE_
      end if

      select case (met_method)
         case (1)
            tausx = ramp*tausx_const
            tausy = ramp*tausy_const

         case (2)

!           Note (KK): old and new meteo data cannot be read at once
!                      since they might come from different files.
!                      Thus only new data is read and new_meteo is set to true.
!                      first call with sst at t_1, later with sst at t_2

#ifdef SLICE_MODEL
            j = jmax/2
#endif
!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)
!$OMP SINGLE

            if (new_meteo) then

               if(calc_met) then

                  if (present(sst_model)) then
! OMP-NOTE: This is an expensive loop, but we cannot thread it as long
!    as exchange_coefficients() and fluxes() pass information through
!    scalars in the meteo module. BJB 2009-09-30.
                     do j=jmin,jmax
                        do i=imin,imax
                           if (az(i,j) .ge. 1) then
                              call exchange_coefficients( &
                                     u10(i,j),v10(i,j),t2(i,j),airp(i,j), &
                                     sst_model(i,j),hum(i,j),hum_method)
                              call fluxes(latc(i,j),u10(i,j),v10(i,j),    &
                                      t2(i,j),tcc(i,j),sst_model(i,j),precip(i,j), &
                                      shf(i,j),tausx(i,j),tausy(i,j),evap(i,j))
                           end if
                        end do
                     end do
                  else
! OMP-NOTE: w needs to be a local (stack) variable to thread this loop.
                     do j=jmin,jmax
                        do i=imin,imax
                           if (az(i,j) .ge. 1) then
! BJB-TODO: Update constants to double.
                              w=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
                              tausx(i,j) = 1.25e-3*1.25*w*U10(i,j)
                              tausy(i,j) = 1.25e-3*1.25*w*V10(i,j)
                           end if
                        end do
                     end do
                  end if

               end if


               call update_2d_halo(airp,airp,az,imin,jmin,imax,jmax,H_TAG)
               call wait_halo(H_TAG)
               call update_2d_halo(tausx,tausx,az,imin,jmin,imax,jmax,H_TAG)
               call wait_halo(H_TAG)
               call update_2d_halo(tausy,tausy,az,imin,jmin,imax,jmax,H_TAG)
               call wait_halo(H_TAG)

               airp_old=>airp_new;airp_new=>airp;airp=>d_airp;d_airp=>airp_old
               tausx_old=>tausx_new;tausx_new=>tausx;tausx=>d_tausx;d_tausx=>tausx_old
               tausy_old=>tausy_new;tausy_new=>tausy;tausy=>d_tausy;d_tausy=>tausy_old
               shf_old=>shf_new;shf_new=>shf;shf=>d_shf;d_shf=>shf_old
               if (calc_met) then
                  tcc_old=>tcc_new;tcc_new=>tcc;tcc=>d_tcc;d_tcc=>tcc_old
               else
                  swr_old=>swr_new;swr_new=>swr;swr=>d_swr;d_swr=>swr_old
               end if
               if (fwf_method .ge. 2) then
                  evap_old=>evap_new;evap_new=>evap;evap=>d_evap;d_evap=>evap_old
               end if
               if (fwf_method.eq.2 .or. fwf_method.eq.3) then
                  precip_old=>precip_new;precip_new=>precip;precip=>d_precip;d_precip=>precip_old
               end if
               if (nudge_sst) then
                  sst_old=>sst_new;sst_new=>sst;sst=>d_sst;d_sst=>sst_old
               end if


               if (.not. first) then
!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
                  do j=jmin-HALO,jmax+HALO
#endif
                     do i=imin-HALO,imax+HALO
                        if (az(i,j) .ne. 0) then
                           d_airp (i,j) = airp_new (i,j) - airp_old (i,j)
                           d_tausx(i,j) = tausx_new(i,j) - tausx_old(i,j)
                           d_tausy(i,j) = tausy_new(i,j) - tausy_old(i,j)
                           d_shf  (i,j) = shf_new  (i,j) - shf_old  (i,j)
                           if (calc_met) then
                              d_tcc(i,j) = tcc_new(i,j) - tcc_old(i,j)
                           else
                              d_swr(i,j) = swr_new(i,j) - swr_old(i,j)
                           end if
                           if (fwf_method .ge. 2) then
                              d_evap(i,j) = evap_new(i,j) - evap_old(i,j)
                           end if
                           if (fwf_method.eq.2 .or. fwf_method.eq.3) then
                              d_precip(i,j) = precip_new(i,j) - precip_old(i,j)
                           end if
                           if (nudge_sst) then
                              d_sst(i,j) = sst_new(i,j) - sst_old(i,j)
                           end if
                        end if
                     end do
#ifndef SLICE_MODEL
                  end do
#endif
!$OMP END DO
!$OMP SINGLE
               end if

            end if


            if (.not. first) then

!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
               do j=jmin-HALO,jmax+HALO
#endif
                  do i=imin-HALO,imax+HALO
                     if (az(i,j) .ne. 0) then
                        airp (i,j) = airp_new (i,j) + d_airp (i,j)*deltm1*t_minus_t2
                        tausx(i,j) = ramp*(tausx_new(i,j) + d_tausx(i,j)*deltm1*t_minus_t2)
                        tausy(i,j) = ramp*(tausy_new(i,j) + d_tausy(i,j)*deltm1*t_minus_t2)
                        shf  (i,j) = shf_new  (i,j) + d_shf  (i,j)*deltm1*t_minus_t2
                        if (calc_met) then
                           tcc(i,j) = tcc_new(i,j) + d_tcc(i,j)*deltm1*t_minus_t2
                           swr(i,j) = short_wave_radiation(yearday,hh,latc(i,j),lonc(i,j),tcc(i,j))
                        else
                           swr(i,j) = swr_new(i,j) + d_swr(i,j)*deltm1*t_minus_t2
                        end if
                        if (fwf_method .ge. 2) then
                           evap(i,j) = evap_new(i,j) + d_evap(i,j)*deltm1*t_minus_t2
                        end if
                        if (fwf_method.eq.2 .or. fwf_method.eq.3) then
                           precip(i,j) = precip_new(i,j) + d_precip(i,j)*deltm1*t_minus_t2
                        end if
                        if (nudge_sst) then
                           sst(i,j) = sst_new(i,j) + d_sst(i,j)*deltm1*t_minus_t2
                        end if
                     end if
                  end do
#ifndef SLICE_MODEL
               end do
#endif
!$OMP END DO
!$OMP SINGLE
#ifdef SLICE_MODEL
               airp (:,j+1) = airp (:,j)
               tausx(:,j+1) = tausx(:,j)
               tausy(:,j+1) = tausy(:,j)
               shf  (:,j+1) = shf  (:,j)
               swr  (:,j+1) = swr  (:,j)
               if (fwf_method .ge. 2) then
                  evap(:,j+1) = evap(:,j)
               end if
               if (fwf_method.eq.2 .or. fwf_method.eq.3) then
                  precip(:,j+1) = precip(:,j)
               end if
               if (nudge_sst) then
                  sst(:,j+1) = sst(:,j)
               end if
#endif

            end if
!$OMP END PARALLEL

         case default

            FATAL 'A non valid meteo method has been specified.'
            stop 'do_meteo'

      end select

   end if
   first = .false.

   call toc(TIM_METEO)
#ifdef DEBUG
     write(debug,*) 'Leaving do_meteo()'
     write(debug,*)
#endif
   return
   end subroutine do_meteo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_meteo() - cleanup the emph{meteo} module.
!
! !INTERFACE:
   subroutine clean_meteo()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine cleans up the \emph{meteo} module.
!
! !REVISION HISTORY:
!  See module for log.
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_meteo() # ',Ncall
#endif

#ifdef DEBUG
     write(debug,*) 'Leaving clean_meteo()'
     write(debug,*)
#endif
   return
   end subroutine clean_meteo
!EOC

!-----------------------------------------------------------------------

   end module meteo

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
