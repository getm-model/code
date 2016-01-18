#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: meteo - provides meteorological forcing for \emph{getm}.
!
! !INTERFACE:
   module getm_ice
!
! !DESCRIPTION:
!
! !SEE ALSO:
!
! !USES:
   use time, only: julianday, secondsofday, timestep
   use domain, only: imin,imax,jmin,jmax,kmax,az
   use domain, only: lonc,latc
   use meteo, only: calc_met
   use meteo, only: shf,swr,albedo,precip,evap,tcc,t2,airp,hum,u10,v10, &
                    tausx,tausy
   use parameters, only: rho_0
   use variables_3d, only: rho,S,T,hn
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
   use ice_winton, only: init_ice_winton, do_ice_winton
   use exceptions
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_getm_ice,do_getm_ice
!
! !PUBLIC DATA MEMBERS:
   integer, public, parameter          :: NO_ICE=0
   integer, public, parameter          :: ICE_THERMODYNAMIC=1
   integer, public, parameter          :: ICE_FROMFILE=2
   integer, public                     :: ice_method=NO_ICE
   integer, public, parameter          :: ICE_FREEZINGPOINT=1
   integer, public, parameter          :: ICE_WINTON=2
   integer, public, parameter          :: ICE_UVIC=3
   integer, public                     :: ice_model=NO_ICE
   character(LEN = PATH_MAX), public   :: ice_file
!  Freezing point ice 'model'
   REALTYPE, public, dimension(:,:), pointer              :: ice_mask=>null()
!  Winton ice model
   REALTYPE, public, dimension(:,:), pointer              :: ice_hi
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_hs
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_ts
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_T1,ice_T2
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_tmelt
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_bmelt
!
! !PUBLIC DATA MEMBERS:
   REALTYPE                            :: hi_thresh=0.2d0
!
! !DEFINED PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_ice - initialise the \emph{meteo} module.
!
! !INTERFACE:
   subroutine init_getm_ice(hotstart)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Basically reads the namelist \emph{meteo} from unit NAMLST.  According to
!  the content of the namelist various variables are allocated and initialised.
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: hotstart
!
! !REVISION HISTORY:
!
!  See log for module.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,rc
   REALTYPE :: ks, alb_sno, alb_ice, pen_ice, opt_dep_ice, opt_ext_ice, &
               opt_ext_snow, t_range_melt, h_lo_lim, kmelt, t_range_dhdt
   namelist /ice/ ice_method,ice_model,ice_file,hi_thresh, &
                  ks, alb_sno, alb_ice, pen_ice, opt_dep_ice, &
                  opt_ext_ice, opt_ext_snow, t_range_melt, h_lo_lim, kmelt, &
                  t_range_dhdt
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_getm_ice'

   ! initialize namelist variables to reasonable defaults
   ks           = 0.31
   alb_sno      = 0.85
   alb_ice      = 0.5826
   pen_ice      = 0.3
   opt_dep_ice  = 0.67
   opt_ext_ice  = 1.5
   opt_ext_snow = 15.0
   t_range_melt = _ONE_
   h_lo_lim     = 0.1
   kmelt        = 240.0
   t_range_dhdt = 0.1

   read(NAMLST,ice)

   select case (ice_method)
      case (NO_ICE)
         LEVEL2 'No ice model included'
         ice_model = NO_ICE
      case (ICE_THERMODYNAMIC)
         !LEVEL2 'Thermodynamic ice model included'
      case (ICE_FROMFILE)
         LEVEL2 'Read in ice mask/thickness from file'
         ice_model = NO_ICE
   end select

   select case (ice_model)
      case (NO_ICE) ! No ice model
         LEVEL2 'No ice model included'
         ice_method = NO_ICE
      case (ICE_FREEZINGPOINT) ! Salinity dependent freezing point
         LEVEL2 'Freezing point ice model'
      case (ICE_WINTON) ! Winton
         LEVEL2 'Winton ice model'
         if (.not. calc_met) then
            call getm_error("init_getm_ice()", "Winton ice model "  // &
                            "requires u10, v10, t2, hum, tcc allocated "// &
                            "which is not the case for calc_met=F."    )
         end if
!        Set ice model parameters
         call  init_ice_winton(ks, alb_sno, alb_ice, pen_ice, &
               opt_dep_ice, opt_ext_ice, opt_ext_snow, t_range_melt, &
               h_lo_lim, kmelt, t_range_dhdt)

!        Allocates memory for the public data members
         allocate(ice_hs(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hs)'
         allocate(ice_T1(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_T1)'
         allocate(ice_T2(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_T2)'
         allocate(ice_tmelt(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_tmelt)'
         allocate(ice_bmelt(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_bmelt)'
         allocate(ice_ts(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_ts)'

         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .ge. 1) then
                  ice_hs(i,j) = _ZERO_
                  ice_T1(i,j) = _ZERO_
                  ice_T2(i,j) = _ZERO_
                  ice_tmelt(i,j) = _ZERO_
                  ice_bmelt(i,j) = _ZERO_
                  ice_ts(i,j) = _ZERO_
               else
                  ice_hs(i,j) = -9999.
                  ice_hi(i,j) = -9999.
                  ice_T1(i,j) = -9999.
                  ice_T2(i,j) = -9999.
                  ice_tmelt(i,j) = -9999.
                  ice_bmelt(i,j) = -9999.
                  ice_ts(i,j) = -9999.
               end if
            end do
         end do
      case default
   end select

   if (ice_method .ne. NO_ICE) then
      allocate(ice_hi(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hi)'
      ice_hi = _ZERO_
      ice_mask => ice_hi
   end if

   return
   end subroutine init_getm_ice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_getm_ice - update the meteo forcing
!
! !INTERFACE:
   subroutine do_getm_ice()
!$ use omp_lib
!
! !DESCRIPTION:
!
!KB   use getm_timers, only: tic, toc, TIM_METEO
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii,jj,n
#if 1
   integer                   :: hum_method=1
   integer                   :: back_radiation_method=1
   integer                   :: fluxes_method=1
#endif
   logical                   :: damp_stress=.false.
   REALTYPE                  :: taucoef=5.0
   REALTYPE                  :: tauh
   REALTYPE                  :: taudamp_avg
   REALTYPE,dimension(I2DFIELD) :: taudamp
!EOP
!-------------------------------------------------------------------------
!BOC
!KB   call tic(TIM_METEO)
   if (ice_method .eq. NO_ICE) return

   if (ice_method .eq. ICE_THERMODYNAMIC) then
   select case (ice_model)
      case (NO_ICE) ! No ice model
      case (ICE_FREEZINGPOINT) ! Salinity dependent freezing point
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .ge. 1 .and. T(i,j,kmax).le.-0.0575*S(i,j,kmax)) then
                  ice_mask(i,j) = _ONE_ 
               else
                  ice_mask(i,j) = _ZERO_ 
               end if
            end do
         end do
      case (2) ! Winton
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .ge. 1) then
                  call do_ice_winton(julianday,secondsofday, &
                         lonc(i,j),latc(i,j), &
                         tcc(i,j),t2(i,j),airp(i,j),hum(i,j), &
                         u10(i,j),v10(i,j), &
                         S(i,j,kmax),rho(i,j,kmax),rho_0,hn(i,j,kmax), &
                         back_radiation_method,hum_method, &
                         fluxes_method,timestep, &
                         T(i,j,kmax),shf(i,j),swr(i,j),precip(i,j), &
                         ice_hs(i,j),ice_hi(i,j),ice_t1(i,j),ice_t2(i,j), &
                         ice_ts(i,j),albedo(i,j),ice_tmelt(i,j),ice_bmelt(i,j))
               else
                  ice_hs(i,j) = -9999.
                  ice_hi(i,j) = -9999.
                  ice_T1(i,j) = -9999.
                  ice_T2(i,j) = -9999.
                  ice_tmelt(i,j) = -9999.
                  ice_bmelt(i,j) = -9999.
                  ice_ts(i,j) = -9999.
              end if
            end do
         end do
      case default
   end select
   end if

   tauh = _HALF_ * hi_thresh
   if (tauh .gt. _ZERO_) then
            call update_2d_halo(ice_hi,ice_hi,az,imin,jmin,imax,jmax,z_TAG)
            call wait_halo(z_TAG)
            taudamp = _ONE_ - _ONE_/(_ONE_ + exp(-taucoef/tauh*(ice_hi-tauh)))
            do j=jmin-HALO+1,jmax+HALO-1
               do i=imin-HALO+1,imax+HALO-1
                  if (az(i,j) .ge. 1) then
                     taudamp_avg = _ZERO_
                     n = 0
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           if (az(ii,jj) .ge. 1) then
                              taudamp_avg = taudamp_avg + taudamp(ii,jj)
                              n = n + 1
                           end if
                        end do
                     end do
                     taudamp_avg = taudamp_avg / n
                     tausx(i,j) = tausx(i,j) * taudamp_avg
                     tausy(i,j) = tausy(i,j) * taudamp_avg
                  end if
               end do
            end do
   else
      tausx = _ZERO_
      tausy = _ZERO_
   end if

!KB    call toc(TIM_METEO)
   return
   end subroutine do_getm_ice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_getm_ice() - cleanup the emph{meteo} module.
!
! !INTERFACE:
   subroutine clean_getm_ice()
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
   return
   end subroutine clean_getm_ice
!EOC

!-----------------------------------------------------------------------

   end module getm_ice

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
