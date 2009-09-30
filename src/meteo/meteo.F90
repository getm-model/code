!$Id: meteo.F90,v 1.21 2009-09-30 11:33:48 bjb Exp $
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
   REALTYPE, public, dimension(:,:), allocatable  :: airp,tausx,tausy,swr,shf
   REALTYPE, public, dimension(:,:), allocatable  :: u10,v10,t2,hum,tcc
   REALTYPE, public, dimension(:,:), allocatable  :: evap,precip
   REALTYPE, public                    :: cd_mom,cd_heat,cd_latent
   REALTYPE, public                    :: cd_precip = _ZERO_
   REALTYPE, public                    :: t_1=-_ONE_,t_2=-_ONE_
   logical, public                     :: new_meteo=.false.
   integer, public                     :: hum_method=-1
!
! !DEFINED PARAMETERS:
   REALTYPE,public,parameter           :: cpa=1008.    !AS that is not exact !
   REALTYPE,public,parameter           :: KELVIN=273.15
   REALTYPE,public,parameter           :: emiss=0.97
   REALTYPE,public,parameter           :: bolz=5.67e-8
!  REALTYPE,public,parameter           :: cpa=1004.67 ! specific heat of dry air- correct
   REALTYPE,public,parameter           :: cpw=4192.   ! specific heat of sea water
   REALTYPE,public,parameter           :: rho_precip = 1000.0
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: meteo.F90,v $
!  Revision 1.21  2009-09-30 11:33:48  bjb
!  Removed tmp output
!
!  Revision 1.20  2009-09-30 11:28:48  bjb
!  OpenMP threading initial implementation
!
!  Revision 1.19  2009-08-18 10:24:46  bjb
!  New getm_timers module
!
!  Revision 1.18  2008-08-29 13:59:12  kb
!  fixed parallel run for constant metforing with convc<>0
!
!  Revision 1.17  2007-06-27 08:39:36  kbk
!  support for fresh water fluxes at the sea surface - Adolf Stips
!
!  Revision 1.16  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.15  2006-03-09 10:53:57  kbk
!  set spinup to -1 when doing hotstart
!
!  Revision 1.14  2006-03-01 13:52:22  kbk
!  renamed method to met_method
!
!  Revision 1.13  2005-04-25 09:25:33  kbk
!  conv --> convc
!
!  Revision 1.12  2005/04/19 15:56:58  kbk
!  added latitude dependent cloud correction factor for long wave rad. - Stips
!
!  Revision 1.11  2005/01/13 09:49:37  kbk
!  wet bulb works, es is global, cleaning - Stips
!
!  Revision 1.10  2004/01/15 11:45:00  kbk
!  meteo point source forcing - taus, swr and shf - implemented
!
!  Revision 1.9  2003/10/01 12:09:13  kbk
!  airp in HALO-zones - need in momentum eqs.
!
!  Revision 1.8  2003/07/01 16:38:34  kbk
!  cleaned code - new methods
!
!  Revision 1.7  2003/06/17 14:53:28  kbk
!  default meteo variables names comply with Adolf Stips suggestion + southpole(3)
!
!  Revision 1.6  2003/05/09 14:28:11  kbk
!  short wave radiation calculated each timestep (not interpolated) - patch from Adolf Stips
!
!  Revision 1.5  2003/04/23 12:05:50  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.4  2003/04/07 15:15:16  kbk
!  merged stable and devel
!
!  Revision 1.3  2003/03/17 15:04:14  gotm
!  Fixed Kondo coefficients - -DWRONG_KONDO can be used
!
!  Revision 1.2  2002/08/16 12:11:06  gotm
!  Fixed parameter order in call to short_wave_radiation()
!
!  Revision 1.1.1.1  2002/05/02 14:01:38  gotm
!  recovering after CVS crash
!
!  Revision 1.8  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.7  2001/09/28 12:32:58  bbh
!  Normalize calculated stresses with 1000. - should be rho_0
!
!  Revision 1.6  2001/08/31 15:41:49  bbh
!  ifdef for CONSTANCE added
!
!  Revision 1.5  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.4  2001/06/04 13:15:12  bbh
!  Further steps towards full implementation of meteorological forcing
!
!  Revision 1.3  2001/05/25 19:14:53  bbh
!  Towards real meteorological forcing
!
!  Revision 1.2  2001/04/24 08:31:18  bbh
!  Initialise all allocated variables with _ZERO_
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: spinup=0,metfmt=2
   REALTYPE                  :: tx= _ZERO_ ,ty= _ZERO_
   REALTYPE                  :: swr_const= _ZERO_ ,shf_const= _ZERO_
   REALTYPE                  :: evap_const= _ZERO_ ,precip_const= _ZERO_
   REALTYPE, dimension(:,:), allocatable :: airp_old,tausx_old,tausy_old
   REALTYPE, dimension(:,:), allocatable :: d_airp,d_tausx,d_tausy
   REALTYPE, dimension(:,:), allocatable :: tcc_old,swr_old,shf_old
   REALTYPE, dimension(:,:), allocatable :: d_tcc,d_swr,d_shf
   REALTYPE, dimension(:,:), allocatable :: evap_old,precip_old
   REALTYPE, dimension(:,:), allocatable :: d_evap,d_precip
!
! !TO DO:
!  A method for stress calculations without knowledge of SST and meteorological
!  variables ($C_d \rho_a W^2) to be used with depth integrated simulations.
!
! !BUGS:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_meteo - initialise the \emph{meteo} module.
!
! !INTERFACE:
   subroutine init_meteo(hotstart)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Basically reads the namelist \emph{meteo} from unit NAMLST.  According to
!  the content of the namelist various variables are allocated and initialised.
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See log for module.
!
! !LOCAL VARIABLES:
   integer                   :: rc
   namelist /meteo/ metforcing,on_grid,calc_met,met_method,fwf_method, &
                    spinup,metfmt,meteo_file, &
                    tx,ty,swr_const,shf_const,evap_const,precip_const, &
                    precip_factor,evap_factor
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

      if (hotstart .and. spinup .gt. 0) then
         LEVEL2 'hotstart --> spinup=-1'
         spinup=-1
      end if
      if ( spinup .gt. 0) then
         LEVEL2 'Forcing will be spun up over ',spinup,' timesteps'
      end if
   end if

! Allocates memory for the public data members

   allocate(airp(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (airp)'
   airp = _ZERO_

   allocate(tausx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (tausx)'
   tausx = _ZERO_

   allocate(tausy(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (tausy)'
   tausy = _ZERO_

   allocate(swr(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (swr)'
   swr = _ZERO_

   allocate(shf(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo: Error allocating memory (shf)'
   shf = _ZERO_

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

   if (metforcing) then

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

      else
      end if

      allocate(airp_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (airp_old)'
      airp_old = _ZERO_

      allocate(tausx_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausx_old)'
      tausx_old = _ZERO_

      allocate(tausy_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tausy_old)'
      tausy_old = _ZERO_

      allocate(d_airp(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_airp)'
      d_airp = _ZERO_

      allocate(d_tausx(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tausx)'
      d_tausx = _ZERO_

      allocate(d_tausy(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tausy)'
      d_tausy = _ZERO_

      allocate(tcc_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (tcc_old)'
      tcc_old = _ZERO_

      allocate(swr_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (swr_old)'
      swr_old = _ZERO_

      allocate(shf_old(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (shf_old)'
      shf_old = _ZERO_

      allocate(d_tcc(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_tcc)'
      d_tcc = _ZERO_

      allocate(d_swr(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_swr)'
      d_swr = _ZERO_

      allocate(d_shf(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo: Error allocating memory (d_shf)'
      d_shf = _ZERO_

      if (fwf_method .ge. 2) then
         allocate(evap_old(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (evap_old)'
         evap_old = _ZERO_

         allocate(d_evap(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_evap)'
         d_evap = _ZERO_
      end if

      if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
         allocate(precip_old(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (precip_old)'
         precip_old = _ZERO_

         allocate(d_precip(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo: Error allocating memory (d_precip)'
         d_precip = _ZERO_
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
   subroutine do_meteo(n,sst)
   use getm_timers, only: tic, toc, TIM_METEO
!$ use omp_lib
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Should be called once every time step to update the meteorological forcing.
!  \emph{do\_meteo()} is called with two arguments - n is the loop number and
!  the sea surface temperature.
!  The SST is only used in the case where fluxes and stresses are
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
!  It is possible to specify a soft start - via the spinup variable -
!  which is used to calculate a ramp (linearly from 0 to one over spinup
!  time steps).
!  To implement an use a different set of formulae for flux calculations
!  should be a matter of only changing the involved subroutines.
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in)                 :: n
   REALTYPE, optional, intent(inout)   :: sst(I2DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See module for log.
!
! !LOCAL VARIABLES:
   integer, save             :: k=0
   integer                   :: i,j
   REALTYPE                  :: ramp,hh,t,t_frac
   REALTYPE                  :: short_wave_radiation
   REALTYPE                  :: uu,cosconv,vv,sinconv
! BJB-TODO: Make sure that 3.14... is defined as *double* precision
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
! BJB-TODO: Change 180. to 180 (integer)
   REALTYPE, parameter       :: deg2rad=pi/180.
   logical,save              :: first=.true.
   logical                   :: have_sst
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_meteo() # ',Ncall
#endif
   call tic(TIM_METEO)

! OMP-NOTE: The vast majority of CPU time is spent in calls to e.g. 
!   exchange_coefficients, fluxes, and especially short_wave_radiation.
!   It is critical to thread these loops, while the rest is not really 
!   critical. In fact, some of the simpler loops have more memory access 
!   than real computations, and should as such not be threaded.
!     BJB 2009-09-30.

!$OMP PARALLEL DEFAULT(SHARED)                                           &
!$OMP     PRIVATE(i,j,ramp,hh,t,t_frac)                                  &
!$OMP     PRIVATE(uu,cosconv,vv,sinconv,have_sst)


   if (metforcing) then

      t = n*timestep

      if(spinup .gt. 0 .and. k .lt. spinup) then
! BJB-TODO: Replace 1.0 with _ONE_ etc in this file.
         ramp = 1.0*k/spinup
!$OMP MASTER
         k = k + 1
!$OMP END MASTER
      else
         ramp = _ONE_
      end if

      select case (met_method)
         case (1)
! BJB-TODO: Why is this called every time step (even after k=spinup)-
!    It should all be constant in time after that(?)

! OMP-NOTE: Memory copy done in serial
!$OMP MASTER
            airp  =  _ZERO_
            tausx = ramp*tx
            tausy = ramp*ty
!$OMP END MASTER
!$OMP BARRIER
!     Rotation of wind stress due to grid convergence
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-1,jmax+1
               do i=imin-1,imax+1
                  if (convc(i,j) .ne. _ZERO_ .and. az(i,j) .gt. 0) then
                     sinconv=sin(-convc(i,j)*deg2rad)
                     cosconv=cos(-convc(i,j)*deg2rad)
                     uu=tausx(i,j)
                     vv=tausy(i,j)
                     tausx(i,j)= uu*cosconv+vv*sinconv
                     tausy(i,j)=-uu*sinconv+vv*cosconv
                  end if
               end do
            end do
!$OMP END DO NOWAIT
!$OMP SINGLE
            swr   = swr_const
            shf   = shf_const
            if (fwf_method .eq. 1) then
               evap = evap_const
               precip = precip_const
            end if
!$OMP END SINGLE
         case (2)
            if(calc_met) then
               have_sst = present(sst)
               if (new_meteo) then
!$OMP MASTER
                  call update_2d_halo(airp,airp,az, &
                                      imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
                  if (.not. first) then
                     tausx_old = tausx
                     tausy_old = tausy
                     tcc_old = tcc
                     shf_old = shf
                     if (fwf_method .ge. 2) then
                        evap_old = evap
                     end if
                     if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                        precip_old = precip
                     end if
                  end if
!$OMP END MASTER
!$OMP BARRIER

                  if (have_sst) then
! BJB-TODO: Check if exchange_coefficients and fluxes are thread safe !$OMP DO SCHEDULE(RUNTIME)
! BJB-TODO: Not tested with sst meteo(!)
!$OMP DO SCHEDULE(RUNTIME)
                     do j=jmin,jmax
                        do i=imin,imax
                           if (az(i,j) .ge. 1) then
                              call exchange_coefficients( &
                                     u10(i,j),v10(i,j),t2(i,j),airp(i,j), &
                                     sst(i,j),hum(i,j),hum_method)
                              call fluxes(latc(i,j),u10(i,j),v10(i,j),    &
                                      t2(i,j),tcc(i,j),sst(i,j),precip(i,j), &
                                      shf(i,j),tausx(i,j),tausy(i,j),evap(i,j))
                           else
                              shf(i,j) = _ZERO_
                              tausx(i,j) = _ZERO_
                              tausy(i,j) = _ZERO_
                              if (fwf_method .ge. 1) then
                                 evap(i,j) = _ZERO_
                                 precip(i,j) = _ZERO_
                              end if
                           end if
                        end do
                     end do
!$OMP END DO
                  else
!$OMP DO SCHEDULE(RUNTIME)
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
!$OMP END DO
                  end if
!$OMP MASTER
                  call update_2d_halo(tausx,tausx,az, &
                                      imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
                  call update_2d_halo(tausy,tausy,az, &
                                      imin,jmin,imax,jmax,H_TAG)
                  call wait_halo(H_TAG)
 !$OMP END MASTER
                  if (.not. first) then
 !$OMP MASTER
                     d_tausx = tausx - tausx_old
                     d_tausy = tausy - tausy_old
                     d_tcc = tcc - tcc_old
                     d_shf = shf - shf_old
                     if (fwf_method .ge. 2) then
                        d_evap = evap - evap_old
                     end if
                     if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                        d_precip = precip - precip_old
                     end if
!$OMP END MASTER
                  end if
               end if
               if (.not. first) then
!$OMP MASTER
                   t_frac = (t-t_1)/(t_2-t_1)
                   shf   = shf_old   + t_frac*d_shf
                   tausx = tausx_old + t_frac*d_tausx
                   tausy = tausy_old + t_frac*d_tausy
                   if (fwf_method .ge. 2) then
                      evap = evap_old + t_frac*d_evap
                   end if
                   if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                      precip = precip_old + t_frac*d_precip
                   end if
!$OMP END MASTER
               end if
! BJB-TODO: Convert constant to full precision:
               hh = secondsofday/3600.
! OMP-NOTE: short_wave_radiation seems thread-safe.
!     Most of the CPU-time is used here, so it is important to thread well.
!$OMP DO SCHEDULE(RUNTIME)
               do j=jmin,jmax
                  do i=imin,imax
                     if (az(i,j) .ge. 1) then
                        swr(i,j) = short_wave_radiation             &
                                (yearday,hh,latc(i,j),lonc(i,j),tcc(i,j))
                     end if
                  end do
               end do
!$OMP END DO
            else
               if (first) then
!$OMP MASTER
                  tausx_old = tausx
                  tausy_old = tausy
                  swr_old = swr
                  shf_old = shf
                  if (fwf_method .ge. 2) then
                     evap_old = evap
                  end if
                  if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                     precip_old = precip
                  end if
!$OMP END MASTER
               end if
               if (new_meteo) then
!$OMP MASTER
                  tausx_old = tausx_old + d_tausx
                  tausy_old = tausy_old + d_tausy
                  swr_old = swr_old + d_swr
                  shf_old = shf_old + d_shf

                  d_tausx = tausx - tausx_old
                  d_tausy = tausy - tausy_old
                  d_swr = swr - swr_old
                  d_shf = shf - shf_old
                  if (fwf_method .ge. 2) then
                     evap_old = evap_old + d_evap
                     d_evap = evap - evap_old
                  end if
                  if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                     precip_old = precip_old + d_precip
                     d_precip = precip - precip_old
                  end if
!$OMP END MASTER
               end if
               if (.not. first) then
                  t_frac = (t-t_1)/(t_2-t_1)
!$OMP MASTER
                  tausx = tausx_old + t_frac*d_tausx
                  tausy = tausy_old + t_frac*d_tausy
                  swr = swr_old + t_frac*d_swr
                  shf = shf_old + t_frac*d_shf
                  if (fwf_method .ge. 2) then
                     evap = evap_old + t_frac*d_evap
                  end if
                  if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
                     precip = precip_old + t_frac*d_precip
                  end if
!$OMP END MASTER
               end if
            end if
         case default
!$OMP MASTER
            FATAL 'A non valid meteo method has been specified.'
            stop 'do_meteo'
!$OMP END MASTER
      end select

   end if

!$OMP END PARALLEL

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
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See module for log.
!
! !LOCAL VARIABLES:
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
