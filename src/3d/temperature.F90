#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  temperature
!
! !INTERFACE:
   module temperature
!
! !DESCRIPTION:
!
! In this module, the temperature equation is processed by
! reading in the namelist {\tt temp} and initialising the temperature field
! (this is done in {\tt init\_temperature}),
! and calculating the advection-diffusion-equation, which includes
! penetrating short-wave radiation as source term (see {\tt do\_temperature}).
!
! !USES:
   use exceptions
   use domain, only: imin,jmin,imax,kmax,jmax,H,az,dry_z
   use variables_3d, only: T,rad,hn,kmin,A,g1,g2
   use meteo, only: metforcing,met_method,nudge_sst,sst,sst_const
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_temperature, do_temperature, init_temperature_field
!
! !PRIVATE DATA MEMBERS:
   integer                   :: temp_method=1,temp_format=2
   character(len=PATH_MAX)   :: temp_file="t_and_s.nc"
   integer                   :: temp_field_no=1
   character(len=32)         :: temp_name='temp'
   REALTYPE                  :: temp_const=20.
   integer                   :: temp_adv_split=0
   integer                   :: temp_adv_hor=1
   integer                   :: temp_adv_ver=1
   REALTYPE                  :: avmolt = 1.4d-7
   integer                   :: temp_AH_method=0
   REALTYPE                  :: temp_AH_const=1.4d-7
   REALTYPE                  :: temp_AH_Prt=_TWO_
   REALTYPE                  :: temp_AH_stirr_const=_ONE_
   integer                   :: attenuation_method=0,jerlov=1
   character(len=PATH_MAX)   :: attenuation_file="attenuation.nc"
   REALTYPE                  :: A_const=0.58,g1_const=0.35,g2_const=23.0
   REALTYPE                  :: swr_bot_refl_frac=-_ONE_
   REALTYPE                  :: swr_min_bot_frac=0.01
   REALTYPE                  :: sst_nudging_time=-_ONE_
   integer                   :: temp_check=0
   REALTYPE                  :: min_temp=-2.,max_temp=35.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_temperature - initialisation of temperature
! \label{sec-init-temperature}
!
! !INTERFACE:
   subroutine init_temperature()
!
! !DESCRIPTION:
!
! Here, the temperature equation is initialised. First, the namelist
! {\tt temp} is read from {\tt getm.inp}. Then, depending on the
! {\tt temp\_method}, the temperature field is read from a
! hotstart file ({\tt temp\_method}=0), initialised with a constant value
! ({\tt temp\_method}=1), initialised and interpolated
! with horizontally homogeneous
! temperature from a given temperature profile ({\tt temp\_method}=2),
! or read in and interpolated from a 3D netCDF field ({\tt temp\_method}=3).
! Finally, a number of sanity checks are performed for the chosen
! temperature advection schemes.
!
! !USES:
   use advection, only: J7
   use advection_3d, only: print_adv_settings_3d
   use variables_3d, only: deformC_3d,deformX_3d,deformUV_3d,calc_stirr
   use m2d, only: Am_method,AM_LES
   use les, only: les_mode,LES_TRACER,LES_BOTH
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: k,i,j,n,rc
   integer                   :: status
   namelist /temp/ &
            temp_method,temp_const,temp_file,                 &
            temp_format,temp_name,temp_field_no,              &
            temp_adv_split,temp_adv_hor,temp_adv_ver,         &
            avmolt,temp_AH_method,temp_AH_const,temp_AH_Prt,  &
            temp_AH_stirr_const,                              &
            attenuation_method,attenuation_file,jerlov,       &
            A_const,g1_const,g2_const,                        &
            swr_bot_refl_frac, swr_min_bot_frac,              &
            sst_nudging_time,temp_check,min_temp,max_temp
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_temperature() # ',Ncall
#endif

   LEVEL2 'init_temperature()'
   read(NAMLST,temp)

   if (avmolt .lt. _ZERO_) then
      LEVEL3 'setting avmolt to 0.'
      avmolt = _ZERO_
   else
      LEVEL3 'avmolt = ',real(avmolt)
   end if

   call init_temperature_field()

!  Sanity checks for advection specifications
   LEVEL3 'Advection of temperature'
   if (temp_adv_hor .eq. J7) stop 'init_temperature: J7 not implemented yet'
   call print_adv_settings_3d(temp_adv_split,temp_adv_hor,temp_adv_ver,_ZERO_)

   select case (temp_AH_method)
      case(0)
         LEVEL3 'temp_AH_method=0 -> horizontal diffusion of heat disabled'
      case(1)
         LEVEL3 'temp_AH_method=1 -> Using constant horizontal diffusivity of heat'
         if (temp_AH_const .lt. _ZERO_) then
              call getm_error("init_temperature()", &
                         "Constant horizontal diffusivity of heat <0");
         end if
         LEVEL4 real(temp_AH_const)
      case(2)
         LEVEL3 'temp_AH_method=2 -> using LES parameterisation'
         LEVEL4 'Turbulent Prandtl number: ',real(temp_AH_Prt)
         deformC_3d =.true.
         deformX_3d =.true.
         deformUV_3d=.true.
         if (Am_method .eq. AM_LES) then
            les_mode = LES_BOTH
         else
            les_mode = LES_TRACER
         end if
      case(3)
         LEVEL3 'temp_AH_method=3 -> SGS stirring parameterisation'
         if (temp_AH_stirr_const .lt. _ZERO_) then
              call getm_error("init_temperature()", &
                         "temp_AH_stirr_const <0");
         end if
         LEVEL4 'stirring constant: ',real(temp_AH_stirr_const)
         deformC_3d =.true.
         deformX_3d =.true.
         deformUV_3d=.true.
         calc_stirr=.true.
      case default
         call getm_error("init_temperature()", &
                         "A non valid temp_AH_method has been chosen");
   end select
   select case (attenuation_method)
      case (0)
         LEVEL3 'setting attenuation coefficients to constant values:'
         select case (jerlov)
            case (1)
               LEVEL4 'Jerlov Type I (=case I): A=0.58, g1=0.35 , g2=23.0'
               A_const=0.58;g1_const=0.35;g2_const=23.0
            case (2)
               LEVEL4 'Jerlov Type 1 (=case I): A=0.68, g1=1.2 , g2=28.0'
               A_const=0.68;g1_const=1.20;g2_const=28.0
            case (3)
               LEVEL4 'Jerlov Type IA (=case IA): A=0.62, g1=0.6 , g2=20.0'
               A_const=0.62;g1_const=0.60;g2_const=20.0
            case (4)
               LEVEL4 'Jerlov Type IB (=case 1B): A=0.67, g1=1.0 , g2=17.0'
               A_const=0.67;g1_const=1.00;g2_const=17.0
            case (5)
               LEVEL4 'Jerlov Type II (=case II): A=0.77, g1=1.5 , g2=14.0'
               A_const=0.77;g1_const=1.50;g2_const=14.0
            case (6)
               LEVEL4 'Jerlov Type III (=case II, coastal): A=0.78, g1=1.40 , g2=7.9'
               A_const=0.78;g1_const=1.40;g2_const=7.9
            case default
               LEVEL4 'User specified:'
               LEVEL4 ' A=  ',A_const
               LEVEL4 ' g1= ',g1_const
               LEVEL4 ' g2= ',g2_const
         end select
         A=A_const
         g1=g1_const
         g2=g2_const
      case (1)
         LEVEL3 'reading attenuation coefficients from:'
         LEVEL4 trim(attenuation_file)
         LEVEL1 'WARNING: reading routine not coded yet'
         LEVEL1 'WARNING: setting to jerlov=1'
         A_const=0.58;g1_const=0.35;g2_const=23.0
         A=A_const
         g1=g1_const
         g2=g2_const
      case default
   end select

   if (swr_bot_refl_frac .gt. _ZERO_) then
      LEVEL3 "reflection of short wave radiation from the bottom:"
      LEVEL4 "swr_bot_refl_frac=",swr_bot_refl_frac
      LEVEL4 "swr_min_bot_frac= ",swr_min_bot_frac
   end if

   if (metforcing .and. sst_nudging_time.gt._ZERO_) then
      nudge_sst = .True.
      LEVEL3 'nudging of SST enabled'
      LEVEL4 'sst_nudging_time=',real(sst_nudging_time)
      allocate(sst(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_temperature: Error allocating memory (sst)'
      select case(met_method)
         case (1)
            if (sst_const .gt. _ZERO_) then
               LEVEL4 'constant sst=',real(sst_const)
               sst = sst_const
            else
               call getm_error("init_temperature()", &
                               "non-positive sst_const")
            end if
         case (2)
            LEVEL4 'sst read from meteo file'
      end select
   end if

   LEVEL3 'temp_check=',temp_check
   if (temp_check .ne. 0) then
      LEVEL4 'doing sanity check on temperature'
      LEVEL4 'min_temp=',min_temp
      LEVEL4 'max_temp=',max_temp
      if (temp_check .gt. 0) then
         LEVEL4 'out-of-bound values result in termination of program'
      end if
      if (temp_check .lt. 0) then
         LEVEL4 'out-of-bound values result in warnings only'
      end if

      call check_3d_fields(imin,jmin,imax,jmax,kmin,kmax,az, &
                           T,min_temp,max_temp,status)
      if (status .gt. 0) then
         if (temp_check .gt. 0) then
            call getm_error("do_temperature()", &
                            "out-of-bound values encountered")
         end if
         if (temp_check .lt. 0) then
            LEVEL1 'do_temperature(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_temperature()'
   write(debug,*)
#endif
   return
   end subroutine init_temperature
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_temperature_field - initialisation of temperature field
! \label{sec-init-temperature-field}
!
! !INTERFACE:
   subroutine init_temperature_field()
!
! !DESCRIPTION:
! Initialise the temperature field as specified with temp\_method
! and exchange the HALO zones
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
   integer                   :: k,i,j,n
   integer, parameter        :: nmax=10000
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: status
!EOP
!-------------------------------------------------------------------------
!BOC

#ifdef NS_FRESHWATER_LENSE_TEST
temp_field_no=1
#endif
#ifdef MED_15X15MINS_TEST
temp_field_no=1
#endif
   select case (temp_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value'
         forall(i=imin:imax,j=jmin:jmax, az(i,j) .ne. 0) &
                T(i,j,:) = temp_const
      case(2)
         LEVEL3 'using profile'
         call read_profile(temp_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,kmax, &
                           az,H,hn,T)
      case(3)
         LEVEL3 'interpolating from 3D field'
         call get_3d_field(temp_file,temp_name,temp_field_no,.true.,T)
      case default
         FATAL 'Not valid temp_method specified'
         stop 'init_temperature'
   end select

   call update_3d_halo(T,T,az,imin,jmin,imax,jmax,kmax,D_TAG)
   call wait_halo(D_TAG)
   call mirror_bdy_3d(T,D_TAG)

#ifdef DEBUG
   write(debug,*) 'Leaving init_temperature_field()'
   write(debug,*)
#endif
   return
   end subroutine init_temperature_field
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_temperature - temperature equation \label{sec-do-temperature}
!
! !INTERFACE:
   subroutine do_temperature(n)
!
! !DESCRIPTION:
!
! Here, one time step for the temperature equation is performed.
! First, preparations for the call to the advection schemes are
! made, i.e.\ calculating the necessary metric coefficients.
! After the call to the advection schemes, which actually perform
! the advection (and horizontal diffusion) step as an operational
! split step, the solar radiation at the interfaces ({\tt rad(k)})
! is calculated
! from given surface radiation ({\tt swr\_loc})
! by means of a double exponential
! approach, see equation (\ref{Light}) on page \pageref{Light}).
! An option to reflect part of the short wave radiation that reaches the
! bottom has been implemented. In very shallow waters - or with very clear
! waters - a significant part of the incoming radiation will reach the
! bottom. Setting swr\_bot\_refl\_frac to a value between 0 and 1 will
! reflect this fraction of what ever the value of SWR is at the bottom.
! The default value of swr\_bot\_refl\_frac is 0.
! The reflection is only done if the ratio between the surface and bottom
! values of SWR is greater than swr\_min\_bot\_frac (default 0.01).
! Furthermore, the surface heat flux {\tt sfl\_loc} is given a
! value.
! The sea surface temperature is limited by the freezing point
! temperature (as a most primitive sea ice model). The next
! step is to set up the tri-diagonal matrix for calculating the
! new temperature by means of a semi-implicit central scheme for the
! vertical diffusion. Source terms which appear on the right hand sides
! are due to the divergence of the solar radiation at the interfaces.
! The subroutine is completed by solving the tri-diagonal linear
! equation by means of a tri-diagonal solver.
!
! !USES:
   use advection_3d, only: do_advection_3d
   use variables_3d, only: dt,cnpar,hn,ho,nuh,uu,vv,ww,hun,hvn,S
   use domain,       only: imin,imax,jmin,jmax,kmax,az
   use meteo,        only: swr,shf
   use parameters,   only: rho_0,cp
   use getm_timers, only: tic,toc,TIM_TEMP,TIM_TEMPH,TIM_MIXANALYSIS
   use variables_3d, only: do_numerical_analyses_3d
   use variables_3d, only: nummix_T,nummix_T_old,nummix_T_int
   use variables_3d, only: phymix_T,phymix_T_int
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: n
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,rc
   REALTYPE                  :: T2(I3DFIELD)
! OMP-NOTE: The pointer declarations is to allow each omp thread to
!   have its own work storage (over a vertical).
   REALTYPE, POINTER         :: Res(:)
   REALTYPE, POINTER         :: auxn(:),auxo(:)
   REALTYPE, POINTER         :: a1(:),a2(:),a3(:),a4(:)
   REALTYPE, POINTER         :: rad1d(:)
   REALTYPE                  :: zz,swr_loc,shf_loc
   REALTYPE                  :: swr_refl
   REALTYPE                  :: rho_0_cpi
   integer                   :: status
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_temperature() # ',Ncall
#endif
   call tic(TIM_TEMP)
   rho_0_cpi = _ONE_/(rho_0*cp)

#ifdef _NUMERICAL_ANALYSES_OLD_
   if (do_numerical_analyses_3d) then
      call toc(TIM_TEMP)
      call tic(TIM_MIXANALYSIS)
! OMP-note: The following array-based line could be implemented
!    with OMP as a WORKSHARE construct. However, it would require a dedicated
!    PARALLEL region (as the various advection schemes have their own regions),
!    so the overhead of the contruct would be rather large.
      T2 = T**2
      call do_advection_3d(dt,T2,uu,vv,ww,hun,hvn,ho,hn,                           &
                           temp_adv_split,temp_adv_hor,temp_adv_ver,_ZERO_,H_TAG)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_TEMP)
   end if
#endif

   call do_advection_3d(dt,T,uu,vv,ww,hun,hvn,ho,hn,                           &
                        temp_adv_split,temp_adv_hor,temp_adv_ver,_ZERO_,H_TAG, &
                        nvd=nummix_T)

#ifdef _NUMERICAL_ANALYSES_OLD_
   if (do_numerical_analyses_3d) then
      call toc(TIM_TEMP)
      call tic(TIM_MIXANALYSIS)
      call numerical_mixing(T2,T,nummix_T_old,nummix_T_int)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_TEMP)
   end if
#endif

   if (temp_AH_method .gt. 0) then
!     T is not halo updated after advection
      call tic(TIM_TEMPH)
      call update_3d_halo(T,T,az,imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)
      call toc(TIM_TEMPH)

      call tracer_diffusion(T,hn,temp_AH_method,temp_AH_const,temp_AH_Prt,temp_AH_stirr_const, &
                            phymix_T)
   end if

   if (do_numerical_analyses_3d) then
      call toc(TIM_TEMP)
      call tic(TIM_MIXANALYSIS)
      call physical_mixing(T,avmolt,phymix_T,phymix_T_int,temp_AH_method)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_TEMP)
   end if

! OMP-NOTE: Pointer definitions and allocation so that each thread can
!           get its own work memory.
!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,rc, zz,swr_loc,shf_loc)                         &
!$OMP    PRIVATE(Res,auxn,auxo,a1,a2,a3,a4,rad1d)

! Each thread allocates its own HEAP storage:
   allocate(Res(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (Res)'
   allocate(auxn(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (auxn)'
   allocate(auxo(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (auxo)'
   allocate(a1(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (a1)'
   allocate(a2(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (a2)'
   allocate(a3(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (a3)'
   allocate(a4(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (auxo)'
   allocate(rad1d(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_temperature: Error allocating memory (rad1d)'

! Note: We do not need to initialize these work arrays.
!   Tested BJB 2009-09-25.

!  Solar radiation and vertical diffusion of temperature

!  Solar radiation
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            swr_loc=swr(i,j)
            rad(i,j,kmax)=swr_loc
            zz = _ZERO_
            do k=kmax-1,0,-1
               zz=zz+hn(i,j,k+1)
               rad(i,j,k)=swr_loc &
                      *(A(i,j)*exp(-zz/g1(i,j))+(1-A(i,j))*exp(-zz/g2(i,j)))
            end do
         end if
      end do
   end do
!$OMP END DO

! OMP-NOTE: This needs local per-thread storage to thread:
!  rad1d, auxo, auxn, a1, a2, a4, a3, Res
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then

            shf_loc=shf(i,j)

            if (T(i,j,kmax).le.-0.0575*S(i,j,kmax)) then  ! use most primitive
                                                          ! sea ice model ...
               shf_loc=max(_ZERO_,shf_loc)
            end if

            do k=0,kmax
               rad1d(k)=rad(i,j,k)
            end do
!           allow for reflection of SWR from bottom.
            if (swr_bot_refl_frac .gt. _ZERO_ .and. &
                rad1d(0)/rad1d(kmax) .gt. swr_min_bot_frac) then
               swr_refl=rad1d(0)*swr_bot_refl_frac
#if 0
            if (D(i,j) .lt. 0.5) then
               STDERR 'SWR ',i,j,swr_loc,swr_refl
            end if
#endif
               rad1d(0)=rad1d(0)-swr_refl
               zz = _ZERO_
               do k=1,kmax
                  zz=zz+hn(i,j,k)
                  rad1d(k)=rad(i,j,k) - swr_refl &
                         *(A(i,j)*exp(-zz/g1(i,j))+(1-A(i,j))*exp(-zz/g2(i,j)))
               end do
            end if
            do k=0,kmax
               rad1d(k)=rad1d(k)*rho_0_cpi                ! note this
            end do

            if (kmax.gt.1) then
!     Auxilury terms, old and new time level,
               do k=1,kmax-1
                  auxo(k)=2*(1-cnpar)*dt*(nuh(i,j,k)+avmolt)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=2*   cnpar *dt*(nuh(i,j,k)+avmolt)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
               end do

!        Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)
               a2(k)=hn(i,j,k)+auxn(k-1)
               a4(k)=T(i,j,k)*(hn(i,j,k)-auxo(k-1))+T(i,j,k-1)*auxo(k-1)  &
                     +dry_z(i,j)*dt*(rad1d(k)+shf_loc*rho_0_cpi-rad1d(k-1))
               if (nudge_sst) then
                  a4(kmax) = a4(kmax) - dry_z(i,j)*dt*hn(i,j,kmax)*(T(i,j,kmax)-sst(i,j))/sst_nudging_time
               end if

!        Matrix elements for inner layers
               do k=2,kmax-1
                  a3(k)=-auxn(k  )
                  a1(k)=-auxn(k-1)
                  a2(k)=hn(i,j,k)+auxn(k)+auxn(k-1)
                  a4(k)=T(i,j,k+1)*auxo(k)                          &
                       +T(i,j,k  )*(hn(i,j,k)-auxo(k)-auxo(k-1))    &
                       +T(i,j,k-1)*auxo(k-1)                        &
                       +dry_z(i,j)*dt*(rad1d(k)-rad1d(k-1))
               end do

!        Matrix elements for bottom layer
               k=1
               a3(k)=-auxn(k  )
               a2(k)=hn(i,j,k)+auxn(k)
               a4(k)=T(i,j,k+1)*auxo(k)                           &
                    +T(i,j,k  )*(hn(i,j,k)-auxo(k))               &
                    +dry_z(i,j)*dt*(rad1d(k)-rad1d(k-1))

               call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)

               do k=1,kmax
                  T(i,j,k)=Res(k)
               end do

            end if
         end if
      end do
   end do
!$OMP END DO

! Each thread must deallocate its own HEAP storage:
   deallocate(Res,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (Res)'
   deallocate(auxn,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (auxn)'
   deallocate(auxo,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (auxo)'
   deallocate(a1,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (a1)'
   deallocate(a2,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (a2)'
   deallocate(a3,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (a3)'
   deallocate(a4,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (a4)'
   deallocate(rad1d,stat=rc)
   if (rc /= 0) stop 'do_temperature: Error deallocating memory (rad1d)'


!$OMP END PARALLEL

   if (temp_check .ne. 0 .and. mod(n,abs(temp_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,kmin,kmax,az, &
                           T,min_temp,max_temp,status)
      if (status .gt. 0) then
         if (temp_check .gt. 0) then
            call getm_error("do_temperature()", &
                            "out-of-bound values encountered")
         end if
         if (temp_check .lt. 0) then
            LEVEL1 'do_temperature(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if


   call toc(TIM_TEMP)
#ifdef DEBUG
   write(debug,*) 'Leaving do_temperature()'
   write(debug,*)
#endif
   return
   end subroutine do_temperature
!EOC

!-----------------------------------------------------------------------

   end module temperature

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
