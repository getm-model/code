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
   use domain, only: imin,imax,jmin,jmax,kmax,az
   use meteo, only: albedo
   use variables_3d, only: S,T
   use ice_winton, only: do_ice_winton
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
   REALTYPE, public, dimension(:,:), allocatable, target  :: ice_T1,ice_T2
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
   integer                   :: rc
   namelist /ice/ ice_method
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_getm_ice'
   read(NAMLST,ice)

   select case (ice_method)
         case (0) ! No ice model
         case (1) ! Salinity dependent freezing point
            allocate(ice_mask(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_mask)'
            ice_mask = _ZERO_
         case (2) ! Winton
!           Allocates memory for the public data members
            allocate(ice_hs(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hs)'
            ice_hs = _ZERO_
            allocate(ice_hi(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_hi)'
            ice_hi = _ZERO_
            allocate(ice_T1(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_T1)'
            ice_T1 = _ZERO_
            allocate(ice_T2(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_getm_ice: Error allocating memory (ice_T2)'
            ice_T2 = _ZERO_
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
#if 0
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .ge. 1) then
                  call do_ice_winton(ice_hs(i,j),ice_hi(i,j),ice_T1(i,j),ice_T2(i,j))
               end if
            end do
         end do
#endif
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
