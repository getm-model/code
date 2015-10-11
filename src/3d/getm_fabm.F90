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
   use domain,only: H
!KB   use get_field, only: get_2d_field,get_3d_field
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d,only: fabm_pel,fabm_ben,fabm_diag,fabm_diag_hz
   use variables_3d, only: nuh,T,S,rho,a,g1,g2,taubmax_3d
   use variables_3d, only: do_numerical_analyses_3d
   use advection_3d, only: print_adv_settings_3d,do_advection_3d
   use variables_2d, only: D,fwf_int
   use meteo, only: swr,wind,evap,precip,tcc
   use time, only: month,yearday,secondsofday,timestr
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
   use exceptions
! JORN_FABM
   use gotm_fabm, only: init_gotm_fabm,set_env_gotm_fabm,do_gotm_fabm
   use gotm_fabm, only: gotm_fabm_calc=>fabm_calc, model, cc_col=>cc, cc_diag_col=>cc_diag, cc_diag_hz_col=>cc_diag_hz, cc_transport

   use fabm, only: type_horizontal_variable_id, fabm_is_variable_used
   use fabm_types,only: output_instantaneous, output_none
   use fabm_standard_variables, only: standard_variables

   IMPLICIT NONE

!
! !PUBLIC DATA MEMBERS:
   public init_getm_fabm, postinit_getm_fabm, do_getm_fabm, model, output_none
   public init_getm_fabm_fields
   integer, public :: fabm_init_method=0
   character(len=PATH_MAX)   :: fabm_init_file
   integer                   :: fabm_init_format, fabm_field_no
   logical, public :: fabm_calc
   REALTYPE,dimension(:,:,:,:),allocatable,target,public :: phymix_fabm_pel,nummix_fabm_pel
!
! !PRIVATE DATA MEMBERS:
   type t_pa3d
      REALTYPE,dimension(:,:,:),pointer :: p3d
   end type t_pa3d
   type(t_pa3d),dimension(:),allocatable :: pa_phymix,pa_nummix
#ifndef _POINTER_REMAP_
   REALTYPE,dimension(:,:,:),allocatable,target :: work3d
#endif
   integer         :: fabm_adv_split=0
   integer         :: fabm_adv_hor=1
   integer         :: fabm_adv_ver=1
   integer         :: fabm_AH_method=0
   REALTYPE        :: fabm_AH_const=1.4d-7
   REALTYPE        :: fabm_AH_Prt=_TWO_
   REALTYPE        :: fabm_AH_stirr_const=_ONE_
   type (type_horizontal_variable_id) :: id_bottom_depth_below_geoid,id_bottom_depth

   type type_input_variable
      integer                              :: ncid  = -1
      integer                              :: varid = -1
      class (type_input_variable), pointer :: next => null()
   end type

   type,extends(type_input_variable) :: type_horizontal_input_variable
      REALTYPE, allocatable, dimension(:,:) :: data
      type (type_horizontal_variable_id)    :: id
   end type

   class (type_input_variable), pointer, save :: first_input_variable => null()

   integer         :: old_month=-1
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

interface
   subroutine tracer_diffusion(f,hn,AH_method,AH_const,AH_Prt,AH_stirr_const, &
                               phymix)
      use domain, only: imin,imax,jmin,jmax,kmax
      IMPLICIT NONE
      REALTYPE,intent(in)           :: hn(I3DFIELD)
      integer,intent(in)            :: AH_method
      REALTYPE,intent(in)           :: AH_const,AH_Prt,AH_stirr_const
      REALTYPE,intent(inout)        :: f(I3DFIELD)
      REALTYPE,dimension(:,:,:),pointer,intent(out),optional :: phymix
   end subroutine tracer_diffusion

   subroutine inquire_file(fn,ncid,varids,varnames)
      character(len=*), intent(in)        :: fn
      integer, intent(inout)              :: ncid
      integer, allocatable, intent(inout) :: varids(:)
      character(len=50), allocatable, intent(out) :: varnames(:)
   end subroutine inquire_file

!KB - only until a proper input_manager has been made
   subroutine get_2d_field_ncdf_by_id(ncid,varid,il,ih,jl,jh,n,field)
      integer, intent(in)                 :: ncid,varid
      integer, intent(in)                 :: il,ih,jl,jh,n
      REALTYPE, intent(out)               :: field(:,:)
   end subroutine get_2d_field_ncdf_by_id

! Temporary interface (should be read from module):
   subroutine get_2d_field(fn,varname,il,ih,jl,jh,break_on_missing,f)
      character(len=*),intent(in)   :: fn,varname
      integer, intent(in)           :: il,ih,jl,jh
      logical, intent(in)           :: break_on_missing
      REALTYPE, intent(out)         :: f(:,:)
   end subroutine get_2d_field
end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_fabm
!
! !INTERFACE:
   subroutine init_getm_fabm(nml_file,hotstart)
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !USES:
   use advection, only: J7
   use variables_3d, only: deformC_3d,deformX_3d,deformUV_3d,calc_stirr
   use m2d, only: Am_method,AM_LES
   use les, only: les_mode,LES_TRACER,LES_BOTH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: nml_file
   logical,intent(in)             :: hotstart
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
   integer, parameter        :: unit_fabm=63
   integer                   :: rc
   integer                   :: i,j,n
   character(len=PATH_MAX)   :: fabm_surface_flux_file=""
   integer                   :: ncid
   integer, allocatable      :: varids(:)
   character(len=50), allocatable :: varnames(:)


   namelist /getm_fabm_nml/ fabm_init_method, &
                           fabm_init_file,fabm_init_format,fabm_field_no, &
                           fabm_surface_flux_file, &
                           fabm_adv_split,fabm_adv_hor,fabm_adv_ver,      &
                           fabm_AH_method,fabm_AH_const,fabm_AH_Prt,      &
                           fabm_AH_stirr_const
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_fabm()'

!  Initialize FABM.
   call init_gotm_fabm(kmax,NAMLST2,'gotm_fabm.nml')

!  Store fabm_calc and model for use by GETM
   fabm_calc = gotm_fabm_calc

   if (fabm_calc) then
      id_bottom_depth_below_geoid = model%get_horizontal_variable_id(standard_variables%bottom_depth_below_geoid)
      id_bottom_depth = model%get_horizontal_variable_id(standard_variables%bottom_depth)

!     Temporary: make sure diagnostic variables store the last value,
!     not their time integral. This will be redundant when time-integrating/averaging
!     is moved from FABM to the physical host.
      do n=1,size(model%diagnostic_variables)
         if (model%diagnostic_variables(n)%output/=output_none) &
            model%diagnostic_variables(n)%output = output_instantaneous
      end do
      do n=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(n)%output/=output_none) &
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
      call print_adv_settings_3d(fabm_adv_split,fabm_adv_hor,fabm_adv_ver,fabm_AH_const)

      select case (fabm_AH_method)
         case(0)
            LEVEL3 'fabm_AH_method=0 -> horizontal diffusion disabled'
         case(1)
            LEVEL3 'fabm_AH_method=1 -> Using constant horizontal diffusivity'
            if (fabm_AH_const .lt. _ZERO_) then
                 call getm_error("init_getm_fabm()", &
                            "Constant horizontal diffusivity <0");
            end if
            LEVEL4 real(fabm_AH_const)
         case(2)
            LEVEL3 'fabm_AH_method=2 -> using LES parameterisation'
            LEVEL4 'Turbulent Prandtl number: ',real(fabm_AH_Prt)
            deformC_3d =.true.
            deformX_3d =.true.
            deformUV_3d=.true.
            if (Am_method .eq. AM_LES) then
               les_mode = LES_BOTH
            else
               les_mode = LES_TRACER
            end if
         case(3)
            LEVEL3 'fabm_AH_method=3 -> SGS stirring parameterisation'
            if (fabm_AH_stirr_const .lt. _ZERO_) then
                 call getm_error("init_getm_fabm()", &
                            "fabm_AH_stirr_const <0");
            end if
            LEVEL4 'stirring constant: ',real(fabm_AH_stirr_const)
            deformC_3d =.true.
            deformX_3d =.true.
            deformUV_3d=.true.
            calc_stirr=.true.
         case default
            call getm_error("init_getm_fabm()", &
                            "A non valid fabm_AH_method has been chosen");
      end select

      allocate(pa_phymix(size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (pa_phymix)'
      allocate(pa_nummix(size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'init_getm_fabm: Error allocating memory (pa_nummix)'
      do n=1,size(model%state_variables)
         pa_phymix(n)%p3d => null()
         pa_nummix(n)%p3d => null()
      end do

!     Here we need to open the NetCDF file with FABM forcing data (if it exists)
!     and loop over all its variables.
!     For now 2D only, so make sure each NetCDF variable is 2D.
!     For each variable, register_horizontal_input_variable should be called (see below).
!     That looks up the FABM variable and also allocates the 2D field that will hold the input data.

      LEVEL2 'FABM input and forcing ...'
      if (len_trim(fabm_surface_flux_file) .ne. 0) then
         LEVEL3 'reading surface fluxes from:'
         LEVEL4 trim(fabm_surface_flux_file)
         call inquire_file(fabm_surface_flux_file,ncid,varids,varnames)
         do n=1,size(varids)
            if ( varids(n) .ne. -1) then
!              remeber surface_flux model in fabm.yaml
               LEVEL4  'inquiring: ',trim(varnames(n))//'_flux'
               call register_horizontal_input_variable(trim(varnames(n))//'_flux',ncid,varids(n))
            end if
         end do
      end if

!     Initialize biogeochemical state variables.
      if (.not. hotstart) then
         call init_getm_fabm_fields()
      end if

   end if

   end subroutine init_getm_fabm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_fabm_fields - initialisation of the fabm fields
! \label{sec-init-getm-fabm-fields}
!
! !INTERFACE:
   subroutine init_getm_fabm_fields()
!
! !DESCRIPTION:
! Initialisation of the getm-fabm fields as specified by fabm\_init\_method
! and exchange of the HALO zones
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
!EOP
!-------------------------------------------------------------------------
!BOC

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
            stop 'init_getm_fabm_fields()'
      end select

!     Update halos with biogeochemical variable values (distribute initial values).
      do n=1,size(model%state_variables)
         call update_3d_halo(fabm_pel(:,:,:,n),fabm_pel(:,:,:,n),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      end do

#ifdef DEBUG
   write(debug,*) 'Leaving init_getm_fabm_fields()'
   write(debug,*)
#endif
   return
   end subroutine init_getm_fabm_fields
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: register_horizontal_input_variable
!
! !INTERFACE:
   subroutine register_horizontal_input_variable(name,ncid,varid)
!
! !DESCRIPTION:
!  Registers FABM horizontal fluxes (surface)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(in) :: name
   integer,         intent(in) :: ncid,varid
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
      class (type_horizontal_input_variable), pointer :: variable
!EOP
!-------------------------------------------------------------------------
!BOC
!  Create the input variable and set associated data (FABM id, 
!  NetCDF id, 2D data field).
   allocate(variable)
   variable%id = model%get_horizontal_variable_id(name)
   if (.not.fabm_is_variable_used(variable%id)) then
      LEVEL2 'Prescribed input variable '//trim(name)//' is not used by FABM.'
      stop 'register_horizontal_input_variable: unrecognized variable name'
   end if
   variable%ncid  = ncid
   variable%varid = varid
   allocate(variable%data(I2DFIELD))
   variable%data = _ZERO_

!  Prepend to the list of inout variables.
   variable%next => first_input_variable
   first_input_variable => variable
   end subroutine register_horizontal_input_variable
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
   REALTYPE,dimension(:,:,:),pointer :: p3d
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

      allocate(phymix_fabm_pel(I3DFIELD,size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'postinit_getm_fabm: Error allocating memory (phymix_fabm_pel)'
      phymix_fabm_pel = _ZERO_

      allocate(nummix_fabm_pel(I3DFIELD,size(model%state_variables)),stat=rc)
      if (rc /= 0) stop 'postinit_getm_fabm: Error allocating memory (nummix_fabm_pel)'
      nummix_fabm_pel = _ZERO_

#ifndef _POINTER_REMAP_
      allocate(work3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'postinit_getm_fabm: Error allocating memory (work3d)'
#endif

      do n=1,size(model%state_variables)
#ifdef _POINTER_REMAP_
         p3d => phymix_fabm_pel(:,:,:,n) ; pa_phymix(n)%p3d(imin-HALO:,jmin-HALO:,0:) => p3d
         p3d => nummix_fabm_pel(:,:,:,n) ; pa_nummix(n)%p3d(imin-HALO:,jmin-HALO:,0:) => p3d
#else
         pa_phymix(n)%p3d => work3d
         pa_nummix(n)%p3d => work3d
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
   REALTYPE        :: I_0,taub_nonnorm,cloud
   REALTYPE        :: z(1:kmax)
   REALTYPE,dimension(I2DFIELD) :: work2d
   class (type_input_variable), pointer :: current_input_variable
   integer         :: ncid,varid
   logical         :: some_var_ok=.false.
!EOP
!-----------------------------------------------------------------------
!BOC

   call tic(TIM_GETM_FABM)

!  First update all input fields
   if (month .ne. old_month) then
      old_month = month
      current_input_variable => first_input_variable
      do while (associated(current_input_variable))
         select type (current_input_variable)
            class is (type_horizontal_input_variable)
               ncid  = current_input_variable%ncid
               varid = current_input_variable%varid
               if (ncid .gt. 0 .and. varid .gt. 0) then
                  some_var_ok = .true.
                  call get_2d_field_ncdf_by_id(ncid,varid,ilg,ihg,jlg,jhg,month, &
                                               current_input_variable%data(ill:ihl,jll:jhl))
               end if
         end select
         current_input_variable => current_input_variable%next
      end do
      if (some_var_ok) then
         LEVEL3 timestr,': reading FABM surface fluxes ...',month
      end if
   end if

   do n=1,size(model%state_variables)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (az(i,j) .eq. 1) then
               fabm_pel(i,j,kmax,n) = fabm_pel(i,j,kmax,n)*(_ONE_-fwf_int(i,j)/ho(i,j,kmax))
            end if
         end do
      end do
   end do

!  First we do all the vertical processes
#ifdef SLICE_MODEL
   do j=2,2
#else
   do j=jmin,jmax
#endif
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then

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
            taub_nonnorm = taubmax_3d(i,j)*rho_0

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
                                   bioshade,I_0,cloud,taub_nonnorm,wind(i,j),precip(i,j),evap(i,j), &
                                   z,A(i,j),g1(i,j),g2(i,j),yearday,secondsofday)
            call model%link_horizontal_data(id_bottom_depth_below_geoid,H(i,j))
            call model%link_horizontal_data(id_bottom_depth,D(i,j))

!           Transfer prescribed input variables
            current_input_variable => first_input_variable
            do while (associated(current_input_variable))
               select type (current_input_variable)
               class is (type_horizontal_input_variable)
                  call model%link_horizontal_data(current_input_variable%id,current_input_variable%data(i,j))
               end select
               current_input_variable => current_input_variable%next
            end do

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

         call do_advection_3d(dt,fabm_pel(:,:,:,n),uu,vv,ww,hun,hvn,ho,hn,           &
                              fabm_adv_split,fabm_adv_hor,fabm_adv_ver,_ZERO_,H_TAG, &
                              nvd=pa_nummix(n)%p3d)

#ifndef _POINTER_REMAP_
         if (do_numerical_analyses_3d) then
            nummix_fabm_pel(:,:,:,n) = pa_nummix(n)%p3d
         end if
#endif

         if (fabm_AH_method .gt. 0) then
            call update_3d_halo(fabm_pel(:,:,:,n),fabm_pel(:,:,:,n),az,imin,jmin,imax,jmax,kmax,D_TAG)
            call wait_halo(D_TAG)
            call tracer_diffusion(fabm_pel(:,:,:,n),hn,fabm_AH_method,fabm_AH_const,fabm_AH_Prt,fabm_AH_stirr_const, &
                                  phymix=pa_phymix(n)%p3d)

#ifndef _POINTER_REMAP_
            if (do_numerical_analyses_3d) then
               phymix_fabm_pel(:,:,:,n) = pa_phymix(n)%p3d
            end if
#endif
         end if

         if (do_numerical_analyses_3d) then
            call physical_mixing(fabm_pel(:,:,:,n),_ZERO_,phymix_fabm_pel(:,:,:,n),work2d,fabm_AH_method)
         end if

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
