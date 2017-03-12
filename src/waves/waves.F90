#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: waves
!
! !INTERFACE:
   module waves
!
! !DESCRIPTION:
!
!
! !USES:
   use variables_waves
   use time           , only: write_time_string,timestr
   use parameters     , only: grav => g
   use exceptions
   use halo_zones     , only: update_2d_halo,wait_halo,H_TAG
   use domain         , only: imin,imax,jmin,jmax,kmax,az,H
   use meteo          , only: metforcing,met_method,wind,u10r,v10r
   use getm_timers    , only: tic,toc,TIM_WAVES

   IMPLICIT NONE
   private
!
! !PUBLIC DATA MEMBERS:
   public init_waves,do_waves,uv_waves,uv_waves_3d
   public stokes_drift_3d
   public bottom_friction_waves,wbbl_tauw,wbbl_rdrag

   integer,public,parameter  :: NO_WAVES=0
   integer,public,parameter  :: WAVES_FROMWIND=1
   integer,public,parameter  :: WAVES_FROMFILE=2
   integer,public,parameter  :: WAVES_FROMEXT=3
   integer,public            :: waveforcing_method=NO_WAVES
   integer,public,parameter  :: WAVES_RS=1
   integer,public,parameter  :: WAVES_VF=2
   integer,public,parameter  :: WAVES_NOSTOKES=3
   integer,public            :: waves_method=WAVES_RS
   character(LEN = PATH_MAX),public :: waves_file
   logical,public            :: on_grid=.true.
   integer,public            :: waves_ramp=0
   integer,public,parameter  :: NO_WBBL=0
   integer,public,parameter  :: WBBL_DATA2=1
   integer,public,parameter  :: WBBL_SOULSBY05=2
   integer,public            :: waves_bbl_method=WBBL_DATA2
   logical,public            :: new_waves=.false.
   logical,public            :: new_StokesC=.false.
!  KK-TODO: for computational efficiency this value should be as small as possible
!           (reduces evaluations of hyperbolic functions)
   REALTYPE,public           :: kD_deepthresh
!   REALTYPE,public,parameter :: kD_deepthresh=100*_ONE_ ! errors<1% for less than 85 layers
!   REALTYPE,public,parameter :: kD_deepthresh= 50*_ONE_ ! errors<1% for less than 40 layers
!   REALTYPE,public,parameter :: kD_deepthresh= 25*_ONE_ ! errors<1% for less than 20 layers
!   REALTYPE,public,parameter :: kD_deepthresh= 10*_ONE_ ! errors<1% for less than  8 layers
!
! !PRIVATE DATA MEMBERS:
   REALTYPE                  :: waves_windscalefactor = _ONE_
   REALTYPE                  :: max_depth_windwaves = -_ONE_
   logical                   :: ramp_is_active=.false.
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine stokes_drift_3d(dt,Dveln,hvel,uuEx,vvEx)
         use domain, only: imin,imax,jmin,jmax,kmax
         IMPLICIT NONE
         REALTYPE,intent(in)                     :: dt
         REALTYPE,dimension(I2DFIELD),intent(in) :: Dveln
         REALTYPE,dimension(I3DFIELD),intent(in) :: hvel
         REALTYPE,dimension(I3DFIELD),intent(inout),optional :: uuEx,vvEx
      end subroutine stokes_drift_3d

      subroutine bottom_friction_waves(U1,V1,DU1,DV1,Dvel,u_vel,v_vel,velU,velV,ru,rv,zub,zvb,taubmax)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in)    :: U1,V1,DU1,DV1,Dvel,u_vel,v_vel,velU,velV
         REALTYPE,dimension(E2DFIELD),intent(inout) :: ru,rv,zub,zvb
         REALTYPE,dimension(:,:),pointer,intent(inout),optional :: taubmax
      end subroutine bottom_friction_waves
   end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_waves - initialising WAVES
! \label{sec-init-waves}
!
! !INTERFACE:
   subroutine init_waves(hotstart,runtype)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in) :: hotstart
   integer,intent(in) :: runtype
!
! !DESCRIPTION:
!
! Here, some necessary memory is allocated (in case of the compiler option
! {\tt STATIC}), and information is written to the log-file of
! the simulation.
!
! !LOCAL VARIABLES
   namelist /waves/ waveforcing_method,waves_method,waves_file,on_grid,&
                    waves_windscalefactor,max_depth_windwaves,         &
                    waves_ramp,waves_bbl_method
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_waves() # ',Ncall
#endif

   LEVEL1 'init_waves'

   read(NAMLST,waves)

   select case (waveforcing_method)
      case(NO_WAVES)
         LEVEL2 'no waveforcing'
         return
      case(WAVES_FROMWIND)
         LEVEL2 'waveforcing data derived from wind data'
         if ( .not. metforcing ) call getm_error("init_waves()",       &
                         "metforcing must be active for WAVES_FROMWIND")
         LEVEL3 'waves_windscalefactor = ',real(waves_windscalefactor)
         if ( max_depth_windwaves .lt. _ZERO_) then
            max_depth_windwaves = 99999.0
         else
            LEVEL3 'max_depth_windwaves = ',real(max_depth_windwaves)
         end if
      case(WAVES_FROMFILE)
         LEVEL2 'waveforcing data read from file: ',trim(waves_file)
         LEVEL3 'on_grid = ',on_grid
      case(WAVES_FROMEXT)
         LEVEL2 'waveforcing data written from external'
      case default
         call getm_error("init_waves()", &
                         "no valid waveforcing_method")
   end select

   select case (waves_method)
      case(WAVES_RS)
         LEVEL2 'waves included via Radiation Stress'
      case(WAVES_VF)
         LEVEL2 'waves included via Vortex Force'
      case(WAVES_NOSTOKES)
         LEVEL2 'waves included only via bottom friction'
      case default
         call getm_error("init_waves()", &
                         "no valid waves_method")
   end select

   call init_variables_waves(runtype)

   if (runtype .eq. 1) then
      kD_deepthresh = 10.0d0
   else
      kD_deepthresh = min( max( 10.0d0 , 1.25d0*kmax ) , log(huge(kD_deepthresh)) )
   end if

   waveK = kD_deepthresh / H

   if (waves_ramp .gt. 1) then
      LEVEL2 'waves_ramp=',waves_ramp
      ramp_is_active = .true.
      if (hotstart) then
         LEVEL3 'WARNING: hotstart is .true. AND waves_ramp .gt. 1'
         LEVEL3 'WARNING: .. be sure you know what you are doing ..'
      end if
   end if

   select case (waves_bbl_method)
      case (NO_WBBL)
         LEVEL2 'no wave BBL'
      case (WBBL_DATA2)
         LEVEL2 'wave BBL by DATA2 formula (Soulsby, 1995, 1997)'
      case (WBBL_SOULSBY05)
         LEVEL2 'wave BBL according to Soulsby & Clarke (2005)'
      case default
         call getm_error("init_waves()", &
                         "no valid_waves_bbl_method")
   end select

   return
   end subroutine init_waves
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_waves -
! \label{sec-do-waves}
!
! !INTERFACE:
   subroutine do_waves(n,D)

! !USES:
   use parameters, only: grav=>g
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                      :: n
   REALTYPE,dimension(E2DFIELD),intent(in) :: D
!
! !DESCRIPTION:
!  D should be depth at velocity time stage.
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD) :: waveECm1
   REALTYPE                     :: depth
   REALTYPE,save                :: ramp=_ONE_
   integer                      :: i,j
   REALTYPE,parameter           :: pi=3.1415926535897932384626433832795029d0
   REALTYPE,parameter           :: twopi = _TWO_*pi
   REALTYPE,parameter           :: oneovertwopi=_ONE_/twopi
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_waves() # ',Ncall
#endif

   call tic(TIM_WAVES)

   select case (waveforcing_method)
      case(WAVES_FROMWIND)
         new_waves = .true.
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               if ( az(i,j) .gt. 0 ) then
                  if (wind(i,j) .gt. _ZERO_) then
                     coswavedir(i,j) = u10r(i,j) / wind(i,j)
                     sinwavedir(i,j) = v10r(i,j) / wind(i,j)
                     depth = min( D(i,j) , max_depth_windwaves )
                     waveH(i,j) = wind2waveHeight(waves_windscalefactor*wind(i,j),depth)
                     waveT(i,j) = wind2wavePeriod(waves_windscalefactor*wind(i,j),depth)
                     waveK(i,j) = wavePeriod2waveNumber(waveT(i,j),D(i,j))
                     waveL(i,j) = twopi / waveK(i,j)
                  else
                     coswavedir(i,j) = _ZERO_
                     sinwavedir(i,j) = _ZERO_
                     waveH(i,j) = _ZERO_
                     waveT(i,j) = _ZERO_
                     waveK(i,j) = kD_deepthresh / D(i,j)
                     waveL(i,j) = _ZERO_
                  end if
               end if
            end do
         end do
      case(WAVES_FROMFILE)
         new_waves = .true.
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               if ( az(i,j) .gt. 0 ) then
                  if (waveL(i,j) .gt. _ZERO_) then
                     waveK(i,j) = twopi / waveL(i,j)
                     waveT(i,j) = twopi / sqrt( grav*waveK(i,j)*tanh(waveK(i,j)*D(i,j)) )
                  else
                     waveK(i,j) = kD_deepthresh / D(i,j)
                     waveT(i,j) = _ZERO_
                  end if
               end if
            end do
         end do
      case(WAVES_FROMEXT)
         if (new_waves) then
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  if ( az(i,j) .gt. 0 ) then
                     waveL(i,j) = twopi / waveK(i,j)
                  end if
               end do
            end do
         end if
   end select


   if (new_waves .or. ramp_is_active) then

      if (ramp_is_active) then
         if (n .ge. waves_ramp) then
            ramp = _ONE_
            ramp_is_active = .false.
            STDERR LINE
            call write_time_string()
            LEVEL3 timestr,': finished waves_ramp=',waves_ramp
            STDERR LINE
         else
            ramp = sqrt(_ONE_*n/waves_ramp)
            waveH = ramp * waveH
         end if
      end if

      waveE = grav * (_QUART_*waveH)**2

!     Note (KK): the stokes_drift routines will still be called, but
!                with zeros. [U|V]StokesC[int|adv] read from a restart
!                file can be nonzero within the first 3d time step!
      if (waves_method .ne. WAVES_NOSTOKES) then
!        calculate depth-integrated Stokes drift at T-point
         waveECm1 = waveE * waveK * waveT * oneovertwopi
         UStokesC = coswavedir * waveECm1
         VStokesC = sinwavedir * waveECm1
         new_StokesC = .true.
      end if

   end if

   new_waves = .false.

   call toc(TIM_WAVES)

#ifdef DEBUG
   write(debug,*) 'Leaving do_waves()'
   write(debug,*)
#endif
   return
   end subroutine do_waves
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_waves -
! \label{sec-uv-waves}
!
! !INTERFACE:
   subroutine uv_waves(UEuler,VEuler,UStokes,VStokes,UStokesC,VStokesC,Dvel,DU,DV,UEx,VEx)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)    :: UEuler,VEuler
   REALTYPE,dimension(E2DFIELD),intent(in)    :: UStokes,VStokes
   REALTYPE,dimension(E2DFIELD),intent(in)    :: UStokesC,VStokesC
   REALTYPE,dimension(E2DFIELD),intent(in)    :: Dvel,DU,DV
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: UEx,VEx
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_waves() # ',Ncall
#endif

   call tic(TIM_WAVES)

   select case(waves_method)
      case (WAVES_RS)
         call radiation_stress(Dvel,UEx,VEx)
      case (WAVES_VF)
         call vortex_force(UEuler,VEuler,UStokes,VStokes,UStokesC,VStokesC,DU,DV,UEx,VEx)
   end select

   call toc(TIM_WAVES)

#ifdef DEBUG
   write(debug,*) 'Leaving uv_waves()'
   write(debug,*)
#endif
   return
   end subroutine uv_waves
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_waves_3d -
! \label{sec-uv-waves-3d}
!
! !INTERFACE:
   subroutine uv_waves_3d(uuEuler,vvEuler,Dveln,hvel,hun,hvn,uuEx,vvEx)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(I2DFIELD),intent(in)    :: Dveln
   REALTYPE,dimension(I3DFIELD),intent(in)    :: uuEuler,vvEuler,hvel,hun,hvn
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: uuEx,vvEx
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_waves_3d() # ',Ncall
#endif

   call tic(TIM_WAVES)

   select case(waves_method)
      case (WAVES_RS)
         call radiation_stress_3d(Dveln,hvel,uuEx,vvEx)
      case (WAVES_VF)
         call vortex_force_3d(uuEuler,vvEuler,hun,hvn,uuEx,vvEx)
   end select

   call toc(TIM_WAVES)

#ifdef DEBUG
   write(debug,*) 'Leaving uv_waves_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_waves_3d
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wind2waveHeight - estimates significant wave height from wind
!
! !INTERFACE:
   REALTYPE function wind2waveHeight(wind,depth,fetch)

! !USES:
   use parameters, only: grav => g
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)          :: wind,depth
   REALTYPE,intent(in),optional :: fetch
!
! !DESCRIPTION:
!  Calculates significant wave height (Hm0).
!  If fetch is not provided, unlimited fetch will be assumed.
!  See page 250 in Holthuijsen (2007).
!
! !LOCAL VARIABLES
   REALTYPE           :: depthstar,fetchstar,waveHeightstar
   REALTYPE           :: wind2,windm2,tanhk3dm3,limiter
   REALTYPE,parameter :: waveHeightstar8 = 0.24d0
   REALTYPE,parameter :: k1 = 4.14d-4
   REALTYPE,parameter :: m1 = 0.79d0
   REALTYPE,parameter :: k3 = 0.343d0
   REALTYPE,parameter :: m3 = 1.14d0
   REALTYPE,parameter :: p  = 0.572d0
!
!EOP
!-----------------------------------------------------------------------
!BOC

   wind2 = wind*wind
   windm2 = _ONE_ / wind2

!  dimensionless depth
   depthstar = grav * depth * windm2

   tanhk3dm3 = tanh(k3*depthstar**m3)

   if (present(fetch)) then
!     dimensionless fetch
      fetchstar = grav * fetch * windm2
      limiter = tanh(k1*fetchstar**m1 / tanhk3dm3)
   else
      limiter = _ONE_
   end if

!  dimensionless significant wave height
   waveHeightstar = waveHeightstar8 * (tanhk3dm3*limiter)**p

!  significant wave height
   wind2waveHeight = wind2 * waveHeightstar / grav

   end function wind2waveHeight
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wind2wavePeriod - estimates peak wave period from wind
!
! !INTERFACE:
   REALTYPE function wind2wavePeriod(wind,depth,fetch)

! !USES:
   use parameters, only: grav => g
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)          :: wind,depth
   REALTYPE,intent(in),optional :: fetch
!
! !DESCRIPTION:
!  Calculates peak wave period.
!  If fetch is not provided, unlimited fetch will be assumed.
!  See page 250 in Holthuijsen (2007).
!  The peak wave period can be empirically related to the significant
!  wave period (Holthuijsen Eqs. (4.2.7) and (4.2.9)).
!
! !LOCAL VARIABLES
   REALTYPE           :: depthstar,fetchstar,wavePeriodstar
   REALTYPE           :: windm2,tanhk4dm4,limiter
   REALTYPE,parameter :: wavePeriodstar8 = 7.69d0
   REALTYPE,parameter :: k2 = 2.77d-7
   REALTYPE,parameter :: m2 = 1.45d0
   REALTYPE,parameter :: k4 = 0.10d0
   REALTYPE,parameter :: m4 = 2.01d0
   REALTYPE,parameter :: q  = 0.187d0
!
!EOP
!-----------------------------------------------------------------------
!BOC

   windm2 = _ONE_ / (wind*wind)

!  dimensionless depth
   depthstar = grav * depth * windm2

   tanhk4dm4 = tanh(k4*depthstar**m4)

   if (present(fetch)) then
!     dimensionless fetch
      fetchstar = grav * fetch * windm2
      limiter = tanh(k2*fetchstar**m2 / tanhk4dm4)
   else
      limiter = _ONE_
   end if

!  dimensionless peak wave period
   wavePeriodstar = wavePeriodstar8 * (tanhk4dm4*limiter)**q

!  peak wave period
   wind2wavePeriod = wind * wavePeriodstar / grav

   end function wind2wavePeriod
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wavePeriod2waveNumber - approximates wave number from wave period
!
! !INTERFACE:
   REALTYPE function wavePeriod2waveNumber(period,depth)

! !USES:
   use parameters, only: grav => g
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in) :: period,depth
!
! !DESCRIPTION:
!  x=k*D=kD,y=omega/sqrt(g/D)=omegastar
!  y=sqrt(x*tanh(x)),y(1)=0.8727=omegastar1
!  x'=lg(x),(dx'|dx)=1/(x*ln(10))
!  y'=lg(y),(dy'|dy)=1/(y*ln(10))
!  m'(x)=(dy'|dx')=(dy'|dy)*(dy|dx)*(dx|dx')=x/y*m(x)
!  m(x)=(dy|dx)=0.5*[tanh(x)+x*(1-tanh(x)**2)]/sqrt(x*tanh(x))
!  m(1)=0.677,m'(1)=0.77572=slopestar1
!  y'=lg(y(1))+m'(1)*x' <=> y=y(1)*[x**m'(1)] <=> x=(y/y(1))**(1/m'(1))
!  shallow: y=x       => x<=y(1)**(1/(1  -m'(1)))=0.5449  => y<=0.5449
!  deep   : y=sqrt(x) => x>=y(1)**(1/(0.5-m'(1)))=1.63865 => y>=1.28
!
!  For alternatives see Holthuijsen (2007) page 124
!  (Eckart, 1952 and Fenton, 1988)
!
! !LOCAL VARIABLES
   REALTYPE           :: omega,omegastar,omegastar2,kD
   REALTYPE,parameter :: sqrtgrav_rec = _ONE_/sqrt(grav)
   REALTYPE,parameter :: omegastar1_rec = _ONE_/0.8727d0
   REALTYPE,parameter :: slopestar1_rec = _ONE_/0.77572d0
   REALTYPE,parameter :: one5th = _ONE_/5
   REALTYPE,parameter :: pi=3.1415926535897932384626433832795029d0
!
!EOP
!-----------------------------------------------------------------------
!BOC

   omega = _TWO_ * pi / period ! radian frequency
   omegastar = omega * sqrt(depth) * sqrtgrav_rec ! non-dimensional radian frequency
   omegastar2 = omegastar*omegastar

!!   approximation by Knut
!!   (errors less than 5%)
!   if ( omegastar .gt. 1.28d0 ) then
!!     deep-water approximation
!      kD = omegastar**2
!   else if ( omegastar .lt. 0.5449d0 ) then
!!     shallow-water approximation
!      kD = omegastar
!   else
!!     tangential approximation in loglog-space for full dispersion relation
!      kD = (omegastar1_rec * omegastar) ** slopestar1_rec
!   end if

!  approximation by Soulsby (1997, page 71) (see (18) in Lettmann et al., 2009)
!  (errors less than 1%)
   if ( omegastar .gt. _ONE_ ) then
      kD = omegastar2 * ( _ONE_ + one5th*exp(_TWO_*(_ONE_-omegastar2)) )
   else
      kD = omegastar * ( _ONE_ + one5th*omegastar2 )
   end if

   wavePeriod2waveNumber = kD / depth

   end function wavePeriod2waveNumber
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wbbl_tauw - calculates wave-only bottom stress
!
! !INTERFACE:
   REALTYPE function wbbl_tauw(waveT,waveH,waveK,depth,z0,wbl)

! !USES:
   use parameters, only: avmmol
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)           :: waveT,waveH,waveK,depth,z0
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out),optional :: wbl
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
   REALTYPE           :: kD,Hrms,omegam1,uorb
   logical,save       :: first=.true.
   REALTYPE,save      :: avmmolm1
   REALTYPE,parameter :: sqrthalf=sqrt(_HALF_)
   REALTYPE,parameter :: pi=3.1415926535897932384626433832795029d0
   REALTYPE,parameter :: oneovertwopi=_HALF_/pi
   REALTYPE,parameter :: ar = 0.24d0 ! 0.26d0
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if (first) then
      avmmolm1 = _ONE_ / avmmol
      first = .false.
   end if

   kD = waveK * depth

   if (waveT.gt._ZERO_ .and. kD.lt.kD_deepthresh) then
      Hrms = sqrthalf * waveH
      omegam1 = oneovertwopi * waveT
!     wave orbital velocity amplitude at bottom (peak orbital velocity, ubot in SWAN)
      uorb = _HALF_ * Hrms / ( omegam1*sinh(kD) )

!     KK-TODO: For combined wave-current flow, the decision on
!              turbulent or laminar flow depends on Rew AND Rec.
!              Furthermore, the decision on rough or smooth flow depends
!              on the final taubmax. (Soulsby & Clarke, 2005)
!              However, here we assume always rough turbulent flows.

!     Note (KK): We do not calculate fw alone, because for small
!                uorb this can become infinite.

!     wave friction factor for rough turbulent flow
      !fwr = 1.39d0 * (aorb/z0)**(-0.52d0)
!     wave-only bottom stress
      !tauw = _HALF_ * fw * uorb**2
      wbbl_tauw = _HALF_ * 1.39d0 * (omegam1/z0)**(-0.52d0) * uorb**(2-0.52d0)

!     bbl thickness (Soulsby & Clarke, 2005)
      if (present(wbl)) wbl = max( 12.0d0*z0 , ar*omegam1*sqrt(wbbl_tauw) )

   else
      wbbl_tauw = _ZERO_
      if (present(wbl)) wbl = 12.0d0*z0
   end if

   end function wbbl_tauw
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wbbl_rdrag - calculates mean bottom friction
!
! !INTERFACE:
   REALTYPE function wbbl_rdrag(tauc,tauw,rdragc,vel,depth,wbbl,z0)

! !USES:
   use parameters, only: kappa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in) :: tauc,tauw,rdragc,vel,depth,wbbl,z0
!
! !DESCRIPTION:
!  rough => total stress => includes form drag (Whitehouse, 2000, page 57)
!  also valid for smooth flows:
!  smooth => skin-friction => z0=avmmol/9/taue_vel
!                          => tauc(cds), tauw(fws), wbbl(as) !!!
!
! !LOCAL VARIABLES:
   REALTYPE :: taue_vel,lnT1m1,lnT2,T3m1,A1dT3,A2dT3,sqrtcddT3,cddT32,taum
   REALTYPE,parameter :: DATA2_a1=1.2d0 ! rough (smooth: 9.0d0; Whitehouse, 2000)
   REALTYPE,parameter :: DATA2_n1=3.2d0 ! rough (smooth: 9.0d0; Whitehouse, 2000)
!
!EOP
!-----------------------------------------------------------------------
!BOC

   select case(waves_bbl_method)
      case (WBBL_DATA2) ! DATA2 formula for rough flow (Soulsby, 1995, 1997)
         wbbl_rdrag = (_ONE_ + DATA2_a1*(tauw/(tauc+tauw))*DATA2_n1) * rdragc
      case (WBBL_SOULSBY05) ! Soulsby & Clarke (2005) for rough flow
         taue_vel = ( tauc**2 + tauw**2 ) ** _QUART_
!!        extension by Malarkey & Davies (2012)
!         taue_vel = (tauc**2 + tauw**2 + _TWO_*tauc*tauw*cos(angle))**_QUART_
         lnT1m1 = _ONE_ / log( wbbl / z0 )
         lnT2 = log( depth / wbbl )
         T3m1 = vel / taue_vel
         A1dT3 = _HALF_ * (lnT2-_ONE_) * lnT1m1
         A2dT3 = kappa * lnT1m1
         sqrtcddT3 = sqrt(A1dT3**2 + T3m1*A2dT3) - A1dT3
         cddT32 = sqrtcddT3*sqrtcddT3
         taum = cddT32 * taue_vel*taue_vel
         wbbl_rdrag = taum / vel
   end select

   end function wbbl_rdrag
!EOC
!-----------------------------------------------------------------------

   end module waves

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
