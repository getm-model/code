#ifdef _FABM_
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
   use domain, only: ilg,ihg,jlg,jhg,ill,ihl,jll,jhl
   use domain, only: az,latc,lonc
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d,only: fabm_pel,fabm_ben,fabm_diag,fabm_diag_hz
   use variables_3d, only: nuh,T,S,rho,a,g1,g2,taub
   use variables_3d, only: do_numerical_analyses_3d
   use advection_3d, only: print_adv_settings_3d,do_advection_3d
   use meteo, only: swr,u10,v10,evap,precip,tcc
   use time, only: yearday,secondsofday
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
! JORN_FABM
   use gotm_fabm, only: init_gotm_fabm,set_env_gotm_fabm,do_gotm_fabm
   use gotm_fabm, only: gotm_fabm_calc=>fabm_calc, model, cc_col=>cc, cc_diag_col=>cc_diag, cc_diag_hz_col=>cc_diag_hz, cc_transport
   use fabm_types,only: output_instantaneous

   IMPLICIT NONE

! Temporary interface (should be read from module):
   interface
      subroutine get_2d_field(fn,varname,il,ih,jl,jh,break_on_missing,f)
         character(len=*),intent(in)   :: fn,varname
         integer, intent(in)           :: il,ih,jl,jh
         logical, intent(in)           :: break_on_missing
         REALTYPE, intent(out)         :: f(:,:)
      end subroutine get_2d_field
   end interface
!
! !PUBLIC DATA MEMBERS:
   public init_getm_fabm, postinit_getm_fabm, do_getm_fabm, model
   integer, public :: fabm_init_method=0
   logical, public :: fabm_calc
   REALTYPE,dimension(:,:,:,:),allocatable,target,public :: nummix_fabm_pel
!
! !PRIVATE DATA MEMBERS:
   type t_pa3d
      REALTYPE,dimension(:,:,:),pointer :: p3d
   end type t_pa3d
   type(t_pa3d),dimension(:),allocatable :: pa_nummix
#ifndef _POINTER_REMAP_
   REALTYPE,dimension(:,:,:),allocatable,target :: nummix
#endif
   integer         :: fabm_adv_split=0
   integer         :: fabm_adv_hor=1
   integer         :: fabm_adv_ver=1
   REALTYPE        :: fabm_AH=-_ONE_
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
   use advection, only: J7
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
                           fabm_adv_split,fabm_adv_hor,fabm_adv_ver,fabm_AH
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_fabm()'

!  Initialize FABM.
   call init_gotm_fabm(kmax,NAMLST2,'gotm_fabm.nml')

!  Store fabm_calc and model for use by GETM
   fabm_calc = gotm_fabm_calc

   if (fabm_calc) then
!     Temporary: make sure diagnostic variables store the last value,
!     not their time integral. This will be redundant when time-integrating/averaging
!     is moved from FABM to the physical host.
      do n=1,size(model%diagnostic_variables)
         model%diagnostic_variables(n)%output = output_instantaneous
      end do
      do n=1,size(model%horizontal_diagnostic_variables)
         model%horizontal_diagnostic_variables(n)%output = output_instantaneous
      end do

!     Allocate memory for pelagic state variables.
      allocate(fabm_pel(I3DFIELD,size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (fabm_pel)'
      fabm_pel = _ZERO_

!     Allocate memory for benthic state variables.
      allocate(fabm_ben(I2DFIELD,size(model%bottom_state_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (fabm_ben)'
      fabm_ben = _ZERO_

!     Allocate memory for 3D diagnostic variables.
      allocate(fabm_diag(I3DFIELD,size(model%diagnostic_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (fabm_diag)'
      fabm_diag = _ZERO_

!     Allocate memory for 2D [horizontal-only] diagnostic variables.
      allocate(fabm_diag_hz(I2DFIELD,size(model%horizontal_diagnostic_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (fabm_diag_hz)'
      fabm_diag_hz = _ZERO_

!     Read settings specific to GETM-FABM interaction.
      open(NAMLST2,status='unknown',file=trim(nml_file))
      read(NAMLST2,NML=getm_fabm_nml)
      close(NAMLST2)

!     Show settings specific to GETM-FABM interaction.
      LEVEL2 'Advection of FABM variables'
      if (fabm_adv_hor .eq. J7) stop 'init_getm_fabm: J7 not implemented yet'
      call print_adv_settings_3d(fabm_adv_split,fabm_adv_hor,fabm_adv_ver,fabm_AH)

      allocate(pa_nummix(size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (pa_nummix)'
      do n=1,size(model%state_variables)
         pa_nummix(n)%p3d => null()
      end do

!     Initialize biogeochemical state variables.
      select case (fabm_init_method)
         case(0)
            LEVEL3 'getting initial biogeochemical fields from hotstart file'
         case(1,2)
            LEVEL3 "initial biogeochemical fields from namelists - fabm.nml"
            do j=jmin,jmax
               do i=imin,imax
                  if (az(i,j) .ge. 1 ) then
                     do n=1,size(model%state_variables)
                        fabm_pel(i,j,:,n) = cc_col(:,n)
                     end do
                     do n=1,size(model%bottom_state_variables)
                        fabm_ben(i,j,  n) = cc_col(1,size(model%state_variables)+n)
                     end do
                  end if
               end do
            end do
            if (fabm_init_method .eq. 2) then
               LEVEL3 'now checking initial fields from ',trim(fabm_init_file)
               do n=1,size(model%state_variables)
                  LEVEL4 'inquiring: ',trim(model%state_variables(n)%name)
                  call get_3d_field(fabm_init_file, &
                                 trim(model%state_variables(n)%name), &
                                 fabm_field_no,.false., &
                                 fabm_pel(:,:,:,n))
               end do
               do n=1,size(model%bottom_state_variables)
                  LEVEL4 'inquiring: ',trim(model%bottom_state_variables(n)%name)
                  call get_2d_field(fabm_init_file, &
                                 trim(model%bottom_state_variables(n)%name), &
                                 ilg,ihg,jlg,jhg,.false., &
                                 fabm_ben(ill:ihl,jll:jhl,n))
               end do
            end if
         case default
            FATAL 'Not valid fabm_init_method specified'
            stop 'init_getm_fabm()'
      end select

!     Update halos with biogeochemical variable values (distribute initial values).
      do n=1,size(model%state_variables)
         call update_3d_halo(fabm_pel(:,:,:,n),fabm_pel(:,:,:,n),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      end do

   end if

   end subroutine init_getm_fabm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: postinit_getm_fabm -
!
! !INTERFACE:
   subroutine postinit_getm_fabm()
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'postinit_getm_fabm() # ',Ncall
#endif

   LEVEL1 'postinit_getm_fabm'

   if (do_numerical_analyses_3d) then

      allocate(nummix_fabm_pel(I3DFIELD,size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'postinit_getm_fabm: Error allocating memory (nummix_fabm_pel)'
      nummix_fabm_pel = _ZERO_

#ifndef _POINTER_REMAP_
      allocate(nummix(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'postinit_getm_fabm: Error allocating memory (nummix)'
#endif

      do n=1,size(model%state_variables)
#ifdef _POINTER_REMAP_
         pa_nummix(n)%p3d(imin-HALO:,jmin-HALO:,0:) => nummix_fabm_pel(:,:,:,n)
#else
         pa_nummix(n)%p3d => nummix
#endif
      end do

   end if

   return
   end subroutine postinit_getm_fabm
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
   use getm_timers, only: tic, toc, TIM_GETM_FABM, TIM_ADVECTFABM
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
   REALTYPE        :: wind_speed,I_0,taub_nonnorm,cloud
   REALTYPE        :: z(1:kmax)
!EOP
!-----------------------------------------------------------------------
!BOC

   call tic(TIM_GETM_FABM)

!  First we do all the vertical processes
#ifdef SLICE_MODEL
   do j=2,2
#else
   do j=jmin,jmax
#endif
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then

!           Calculate wind speed from wind vector components.
            if (allocated(u10) .and. allocated(v10)) then
               wind_speed = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
            else
               wind_speed = _ZERO_
            end if

!           Get surface short-wave radiation.
            if (associated(swr)) then
               I_0 = swr(i,j)
            else
               I_0 = _ZERO_
            end if

            if (associated(tcc)) then
               cloud = tcc(i,j)
            else
               cloud = _ZERO_
            end if

!           Calculate depths of cell centers from layer heights.
            z(kmax) = -_HALF_*hn(i,j,kmax)
            do k=kmax-1,1,-1
               z(k) = z(k+1) - _HALF_*(hn(i,j,k+1)+hn(i,j,k))
            end do

!           Calculate actual bottom stress from normalized bottom stress (taub/rho_0)
            taub_nonnorm = taub(i,j)*rho_0

!           Copy current values of biogeochemical variables from full 3D field to columns.
            do n=1,size(model%state_variables)
               cc_col(:,n) = fabm_pel(i,j,:,n)
            end do
            do n=1,size(model%bottom_state_variables)
               cc_col(1,size(model%state_variables)+n) = fabm_ben(i,j,n)
            end do
            do n=1,size(model%diagnostic_variables)
               cc_diag_col(:,n) = fabm_diag(i,j,1:,n)
            end do
            do n=1,size(model%horizontal_diagnostic_variables)
               cc_diag_hz_col(n) = fabm_diag_hz(i,j,n)
            end do

!           Transfer pointers to physical environment variables to FABM.
            call set_env_gotm_fabm(latc(i,j),lonc(i,j),dt,0,0,T(i,j,1:),S(i,j,1:), &
                                   rho(i,j,1:),nuh(i,j,0:),hn(i,j,0:),ww(i,j,0:), &
                                   bioshade,I_0,cloud,taub_nonnorm,wind_speed,precip(i,j),evap(i,j), &
                                   z,A(i,j),g1(i,j),g2(i,j),yearday,secondsofday)

!           Update biogeochemical variable values.
            call do_gotm_fabm(kmax)

!           Copy updated column values of biogeochemical variables to full 3D field.
            do n=1,size(model%state_variables)
               fabm_pel(i,j,:,n) = cc_col(:,n)
            end do
            do n=1,size(model%bottom_state_variables)
               fabm_ben(i,j,n) = cc_col(1,size(model%state_variables)+n)
            end do
            do n=1,size(model%diagnostic_variables)
               fabm_diag(i,j,1:,n) = cc_diag_col(:,n)
            end do
            do n=1,size(model%horizontal_diagnostic_variables)
               fabm_diag_hz(i,j,n) = cc_diag_hz_col(n)
            end do

         end if
      end do
   end do

#ifdef SLICE_MODEL
   do i=imin,imax
      fabm_pel(i,3,:,:)=fabm_pel(i,2,:,:)
      fabm_ben(i,3,:)  =fabm_ben(i,2,:)
   end do
#endif

!  Advect pelagic biogeochemical variables.
   call tic(TIM_ADVECTFABM)
   do n=1,size(model%state_variables)

      if (cc_transport(n)) then
         call update_3d_halo(fabm_pel(:,:,:,n),fabm_pel(:,:,:,n),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)

         call do_advection_3d(dt,fabm_pel(:,:,:,n),uu,vv,ww,hun,hvn,ho,hn,            &
                              fabm_adv_split,fabm_adv_hor,fabm_adv_ver,fabm_AH,H_TAG, &
                              nvd=pa_nummix(n)%p3d)
#ifndef _POINTER_REMAP_
         if (do_numerical_analyses_3d) then
            nummix_fabm_pel(:,:,:,n) = pa_nummix(n)%p3d
         end if
#endif
      end if

   end do
   call toc(TIM_ADVECTFABM)

   call toc(TIM_GETM_FABM)

   end subroutine do_getm_fabm
!EOC

!-----------------------------------------------------------------------

   end module getm_fabm
#endif

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
