#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: getm_fabm()
!
! !INTERFACE:
   module getm_fabm
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn,fadv3d
   use variables_3d, only: nuh,T,S,rho,a,g1,g2,taub
   use advection_3d, only: print_adv_settings_3d,do_advection_3d
   use meteo, only: swr,u10,v10,evap,precip
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
! JORN_FABM
   use gotm_fabm, only: init_gotm_fabm,set_env_gotm_fabm,do_gotm_fabm
   use gotm_fabm, only: fabm_calc, model, cc_col=>cc, cc_diag_col=>cc_diag, cc_diag_hz_col=>cc_diag_hz
   use fabm_types,only: time_treatment_last

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_getm_fabm, do_getm_fabm
   integer, public           :: fabm_init_method=0
!
! !PRIVATE DATA MEMBERS:
   integer         :: fabm_adv_split=0
   integer         :: fabm_hor_adv=1
   integer         :: fabm_ver_adv=1
   REALTYPE        :: fabm_AH=-_ONE_
   REALTYPE, allocatable, dimension(:,:,:,:) :: cc_pel,cc_diag
   REALTYPE, allocatable, dimension(:,:,:)   :: cc_ben,cc_diag_hz
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_fabm
!
! !INTERFACE:
   subroutine init_getm_fabm(nml_file)
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: nml_file
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
   integer, parameter        :: unit_fabm=63
   integer                   :: rc
   integer                   :: i,j,n
   character(len=PATH_MAX)   :: fabm_init_file
   integer                   :: fabm_init_format, fabm_field_no

   namelist /getm_fabm_nml/ fabm_init_method, &
                           fabm_init_file,fabm_init_format,fabm_field_no, &
                           fabm_hor_adv,fabm_ver_adv,fabm_adv_split,fabm_AH
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_fabm()'

!  Initialize FABM.
   call init_gotm_fabm(kmax,NAMLST2,'fabm.nml')

   if (fabm_calc) then
!     Temporary: make sure diagnostic variables store the last value,
!     not their time integral. This will be redundant when time-integrating/averaging
!     is moved from FABM to the physical host.
      do n=1,ubound(model%info%diagnostic_variables,1)
         model%info%diagnostic_variables(n)%time_treatment = time_treatment_last
      end do
      do n=1,ubound(model%info%diagnostic_variables_hz,1)
         model%info%diagnostic_variables_hz(n)%time_treatment = time_treatment_last
      end do

!     Allocate memory for pelagic state variables.
      allocate(cc_pel(ubound(model%info%state_variables,1),I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (cc_pel)'
      cc_pel = _ZERO_

!     Allocate memory for benthic state variables.
      allocate(cc_ben(ubound(model%info%state_variables_ben,1),I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (cc_ben)'
      cc_ben = _ZERO_

!     Allocate memory for 3D diagnostic variables.
      allocate(cc_diag(ubound(model%info%diagnostic_variables,1),I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (cc_diag)'
      cc_diag = _ZERO_

!     Allocate memory for 2D [horizontal-only] diagnostic variables.
      allocate(cc_diag_hz(ubound(model%info%diagnostic_variables_hz,1),I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (cc_diag_hz)'
      cc_diag_hz = _ZERO_

!     Read settings specific to GETM-FABM interaction.
      open(NAMLST2,status='unknown',file=trim(nml_file))
      read(NAMLST2,NML=getm_fabm_nml)
      close(NAMLST2)

!     Show settings specific to GETM-FABM interaction.
      LEVEL2 'Advection of FABM variables'
      call print_adv_settings_3d(fabm_adv_split,fabm_hor_adv,fabm_ver_adv,fabm_AH)

!     Initialize biogeochemical state variables.
      select case (fabm_init_method)
         case(0)
            LEVEL3 'getting initial biogeochemical fields from hotstart file'
         case(1)
            LEVEL3 "initial biogeochemical fields from namelists - fabm.nml"
            do j=jmin,jmax
               do i=imin,imax
                  if (az(i,j) .ge. 1 ) then
                     cc_pel(:,i,j,:) = cc_col(1:ubound(model%info%state_variables,1) ,:)
                     cc_ben(:,i,j)   = cc_col(ubound(model%info%state_variables,1)+1:,1)
                  end if
               end do
            end do
         case(2)
            LEVEL3 'reading initial biogeochemical fields from ',trim(fabm_init_file)
            do n=1,ubound(model%info%state_variables,1)
               LEVEL4 'inquiring ',trim(model%info%state_variables(n)%name)
               call get_field(fabm_init_file,trim(model%info%state_variables(n)%name),fabm_field_no, &
                              cc_pel(n,:,:,:))
            end do
         case default
            FATAL 'Not valid fabm_init_method specified'
            stop 'init_getm_fabm()'
      end select

!     Update halos with biogeochemical variable values (distribute initial values).
      do n=1,ubound(model%info%state_variables,1)
         call update_3d_halo(cc_pel(n,:,:,:),cc_pel(n,:,:,:),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      end do

   end if

   return
   end subroutine init_getm_fabm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_getm_fabm()
!
! !INTERFACE:
   subroutine do_getm_fabm(dt)
!
! !DESCRIPTION:
!
! !USES:
   use getm_timers, only: tic, toc, TIM_GETM_BIO
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer         :: n
   integer         :: i,j,k
   REALTYPE        :: bioshade(1:kmax)
   REALTYPE        :: wind_speed,I_0
   REALTYPE        :: z(1:kmax)
!EOP
!-----------------------------------------------------------------------
!BOC

   call tic(TIM_GETM_BIO)

!  First we do all the vertical processes
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then

!           Calculate wind speed from wind vector components.
            if (allocated(u10) .and. allocated(v10)) then
               wind_speed = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
            else
               wind_speed = _ZERO_
            end if

!           Get surface short-wave radiation.
            if (allocated(swr)) then
               I_0 = swr(i,j)
            else
               I_0 = _ZERO_
            end if

!           Calculate depths of cell centers from layer heights.
            z(kmax) = -_HALF_*hn(i,j,kmax)
            do k=kmax-1,1,-1
               z(k) = z(k+1) - _HALF_*(hn(i,j,k+1)+hn(i,j,k))
            end do

!           Copy current values of biogeochemical variables from full 3D field to columns.
            cc_col(1:ubound(model%info%state_variables,1) ,:) = cc_pel(:,i,j,:)
            cc_col(ubound(model%info%state_variables,1)+1:,1) = cc_ben(:,i,j)
            cc_diag_col    = cc_diag(:,i,j,:)
            cc_diag_hz_col = cc_diag_hz(:,i,j)

!           Transfer pointers to physical environment variables to FABM.
            call set_env_gotm_fabm(dt,0,0,T(i,j,1:),S(i,j,1:), &
                                   rho(i,j,1:),nuh(i,j,0:),hn(i,j,0:),ww(i,j,0:), &
                                   bioshade,I_0,taub(i,j),wind_speed,precip(i,j),evap(i,j), &
                                   z,A(i,j),g1(i,j),g2(i,j))

!           Update biogeochemical variable values.
            call do_gotm_fabm(kmax)

!           Copy updated column values of biogeochemical variables to full 3D field.
            cc_pel    (:,i,j,:) = cc_col(1:ubound(model%info%state_variables,1) ,:)
            cc_ben    (:,i,j)   = cc_col(ubound(model%info%state_variables,1)+1:,1)
            cc_diag   (:,i,j,:) = cc_diag_col
            cc_diag_hz(:,i,j)   = cc_diag_hz_col

         end if
      end do
   end do

!  Advect pelagic biogeochemical variables.
   do n=1,ubound(model%info%state_variables,1)

#if 1
      fadv3d = cc_pel(n,:,:,:)

      call update_3d_halo(fadv3d,fadv3d,az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

!  KK-TODO: fabm_AH_method + include fabm_AH_method=1 into advection

      call do_advection_3d(dt,fadv3d,uu,vv,ww,hun,hvn,ho,hn,                       &
                           fabm_hor_adv,fabm_ver_adv,fabm_adv_split,fabm_AH,H_TAG)

!      if (fabm_AH_method .gt. 1) then
!         call update_3d_halo(fadv3d,fadv3d,az,imin,jmin,imax,jmax,kmax,D_TAG)
!         call wait_halo(D_TAG)
!         call tracer_diffusion(ff,fabm_AH_method,fabm_AH_const,fabm_AH_Prt,fabm_AH_stirr_const)
!      end if

      cc_pel(n,:,:,:) = fadv3d
#else
      call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

!  KK-TODO: fabm_AH_method + include fabm_AH_method=1 into advection

      call do_advection_3d(dt,cc3d(n,:,:,:),uu,vv,ww,hun,hvn,ho,hn,                &
                           fabm_hor_adv,fabm_ver_adv,fabm_adv_split,fabm_AH,H_TAG)

!      if (fabm_AH_method .gt. 1) then
!         call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az,imin,jmin,imax,jmax,kmax,D_TAG)
!         call wait_halo(D_TAG)
!         call tracer_diffusion(cc3d(n,:,:,:),fabm_AH_method,fabm_AH_const,fabm_AH_Prt,fabm_AH_stirr_const)
!      end if
#endif
   end do

   call toc(TIM_GETM_BIO)

   return
   end subroutine do_getm_fabm
!EOC

!-----------------------------------------------------------------------

   end module getm_fabm

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
