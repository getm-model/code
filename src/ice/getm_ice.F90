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
   use meteo, only: shf,swr,albedo,precip,evap,tcc,t2,airp,hum,u10,v10,kelvin
   use parameters, only: rho_0
   use variables_3d, only: rho,rho_0,S,T,hn
   use ice_winton, only: init_ice_winton, do_ice_winton
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_getm_ice,do_getm_ice
!
! !PUBLIC DATA MEMBERS:
   integer, public                     :: ice_method=0
!  Freezing point ice 'model'
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_mask
!  Winton ice model
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_hs,ice_hi
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_ts
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_T1,ice_T2
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_tmelt
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_bmelt
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
   namelist /ice/ ice_method, ks, alb_sno, alb_ice, pen_ice, opt_dep_ice, &
                  opt_ext_ice, opt_ext_snow, t_range_melt, h_lo_lim, kmelt, &
                  t_range_dhdt
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_getm_ice'
   ! initialize namelist variables to reasonable defaults
   ice_method   = 0
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
      case (0) ! No ice model
         LEVEL2 'No ice model included'
      case (1) ! Salinity dependent freezing point
         LEVEL2 'Freezing point ice model'
         allocate(ice_mask(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_mask)'
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .ge. 1) then
                  ice_mask(i,j) = _ZERO_
               else
                  ice_mask(i,j) = -9999.
               end if
            end do
         end do
      case (2) ! Winton
         LEVEL2 'Winton ice model'
!        Set ice model parameters
         call  init_ice_winton(ks, alb_sno, alb_ice, pen_ice, &
               opt_dep_ice, opt_ext_ice, opt_ext_snow, t_range_melt, &
               h_lo_lim, kmelt, t_range_dhdt)
!        Allocates memory for the public data members
         allocate(ice_hs(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hs)'
         allocate(ice_hi(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hi)'
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
                  ice_hi(i,j) = _ZERO_
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
   integer                   :: i,j
#if 1
   integer                   :: hum_method=1
   integer                   :: back_radiation_method=1
   integer                   :: fluxes_method=1
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
!KB   call tic(TIM_METEO)
   select case (ice_method)
      case (0) ! No ice model
      case (1) ! Salinity dependent freezing point
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
                         tcc(i,j),t2(i,j)-kelvin,airp(i,j),hum(i,j), &
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
