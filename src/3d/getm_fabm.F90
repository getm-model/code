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
   use parameters, only: rho_0
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d, only: nuh,T,S,rho,a,g1,g2,taub
   use advection_3d, only: do_advection_3d
   use meteo, only: swr,u10,v10,evap,precip
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
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
   integer         :: fabm_hor_adv=1
   integer         :: fabm_ver_adv=1
   integer         :: fabm_adv_split=1
   REALTYPE        :: fabm_AH=-1.
#ifdef STATIC
   REALTYPE        :: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE        :: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE        :: area_inv(I2DFIELD)
   REALTYPE        :: ff(I3DFIELD)
#else
   REALTYPE, dimension(:,:), allocatable :: delxu,delxv
   REALTYPE, dimension(:,:), allocatable :: delyu,delyv
   REALTYPE, dimension(:,:), allocatable :: area_inv
   REALTYPE, dimension(:,:,:), allocatable :: ff
#endif

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
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
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
      LEVEL2 "settings related to FABM calculations"
      LEVEL3 'fabm_hor_adv=   ',fabm_hor_adv
      LEVEL3 'fabm_ver_adv=   ',fabm_ver_adv
      LEVEL3 'fabm_adv_split= ',fabm_adv_split
      LEVEL3 'fabm_AH=        ',fabm_AH

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

#ifndef STATIC
      allocate(delxu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (delxu)'

      allocate(delxv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (delxv)'

      allocate(delyu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (delyu)'

      allocate(delyv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (delyv)'

      allocate(area_inv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (area_inv)'

      allocate(ff(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (ff)'
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
      delxu=dxu
      delxv=dxv
      delyu=dyu
      delyv=dyv
      area_inv=arcd1
#else
      delxu=dx
      delxv=dx
      delyu=dy
      delyv=dy
      area_inv=ard1
#endif
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
   REALTYPE        :: wind_speed,I_0,taub_nonnorm
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
            
!           Calculate actual bottom stress from normalized bottom stress (taub/rho_0)
            taub_nonnorm = taub(i,j)*rho_0

!           Copy current values of biogeochemical variables from full 3D field to columns.
            cc_col(1:ubound(model%info%state_variables,1) ,:) = cc_pel(:,i,j,:)
            cc_col(ubound(model%info%state_variables,1)+1:,1) = cc_ben(:,i,j)
            cc_diag_col    = cc_diag(:,i,j,:)
            cc_diag_hz_col = cc_diag_hz(:,i,j)

!           Transfer pointers to physical environment variables to FABM.
            call set_env_gotm_fabm(dt,0,0,T(i,j,1:),S(i,j,1:), &
                                   rho(i,j,1:),nuh(i,j,0:),hn(i,j,0:),ww(i,j,0:), &
                                   bioshade,I_0,taub_nonnorm,wind_speed,precip(i,j),evap(i,j), &
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
      ff = cc_pel(n,:,:,:)

      call update_3d_halo(ff,ff,az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

      call do_advection_3d(dt,ff,uu,vv,ww,hun,hvn,ho,hn, &
              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
              fabm_hor_adv,fabm_ver_adv,fabm_adv_split,fabm_AH)

      cc_pel(n,:,:,:) = ff
#else
      call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)

      call do_advection_3d(dt,cc3d(n,:,:,:),uu,vv,ww,hun,hvn,ho,hn, &
              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
              fabm_hor_adv,fabm_ver_adv,fabm_adv_split,fabm_AH)
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
