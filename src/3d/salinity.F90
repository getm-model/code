#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  Salinity
!
! !INTERFACE:
   module salinity
!
! !DESCRIPTION:
!
! In this module, the salinity equation is processed by
! reading in the namelist {\tt salt} and initialising the salinity field
! (this is done in {\tt init\_salinity}),
! and calculating the advection-diffusion-equation (see {\tt do\_salinity}).
!
! !USES:
   use exceptions
   use domain, only: imin,jmin,imax,jmax,kmax,H,az,dry_z
!KB   use get_field, only: get_3d_field
   use variables_2d, only: fwf_int
   use variables_3d, only: S,hn,kmin
   use meteo, only: metforcing,met_method,nudge_sss,sss
   use meteo, only: METEO_CONST,METEO_FROMFILE,METEO_FROMEXT
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
   IMPLICIT NONE

   interface
      subroutine tracer_diffusion(f,hn,AH_method,AH_const,AH_Prt,AH_stirr_const, &
                                  phymix)
         use domain, only: imin,imax,jmin,jmax,kmax
         REALTYPE,intent(in)           :: hn(I3DFIELD)
         integer,intent(in)            :: AH_method
         REALTYPE,intent(in)           :: AH_const,AH_Prt,AH_stirr_const
         REALTYPE,intent(inout)        :: f(I3DFIELD)
         REALTYPE,dimension(:,:,:),pointer,intent(out),optional :: phymix
      end subroutine tracer_diffusion
   end interface
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_salinity, do_salinity, init_salinity_field
!
! !PRIVATE DATA MEMBERS:
   integer                   :: salt_method=1,salt_format=2
   character(len=PATH_MAX)   :: salt_file="t_and_s.nc"
   integer                   :: salt_field_no=1
   character(len=32)         :: salt_name='salt'
   REALTYPE                  :: salt_const=35*_ONE_
   integer                   :: salt_adv_split=0
   integer                   :: salt_adv_hor=1
   integer                   :: salt_adv_ver=1
   REALTYPE                  :: avmols = 1.1d-9
   integer                   :: salt_AH_method=0
   REALTYPE                  :: salt_AH_const=1.1d-9
   REALTYPE                  :: salt_AH_Prt=_TWO_
   REALTYPE                  :: salt_AH_stirr_const=_ONE_
   REALTYPE                  :: sss_nudging_time=-_ONE_
   integer                   :: salt_check=0
   REALTYPE                  :: min_salt=_ZERO_,max_salt=40*_ONE_
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
! !IROUTINE: init_salinity - initialisation of salinity
! \label{sec-init-salinity}
!
! !INTERFACE:
   subroutine init_salinity(hotstart)
!
! !DESCRIPTION:
!
! Here, the salinity equation is initialised. First, the namelist
! {\tt salt} is read from {\tt getm.inp}. Then, depending on the
! {\tt salt\_method}, the salinity field is read from a
! hotstart file ({\tt salt\_method}=0), initialised with a constant value
! ({\tt salt\_method}=1), initialised and interpolated
! with horizontally homogeneous
! salinity from a given salinity profile ({\tt salt\_method}=2),
! or read in and interpolated from a 3D netCDF field ({\tt salt\_method}=3).
! Finally, a number of sanity checks are performed for the chosen
! salinity advection schemes.
!
! Apart from this, there are various options for specific initial
! conditions which are selected by means of compiler options.
!
! !USES:
   use advection, only: J7
   use advection_3d, only: print_adv_settings_3d
   use variables_3d, only: deformC_3d,deformX_3d,deformUV_3d,calc_stirr
   use m2d, only: Am_method,AM_LES
   use les, only: les_mode,LES_TRACER,LES_BOTH
   use meteo, only: sss_const
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in)        :: hotstart
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n,rc
   integer                   :: status
   NAMELIST /salt/                                            &
            salt_method,salt_const,salt_file,                 &
            salt_format,salt_name,salt_field_no,              &
            salt_adv_split,salt_adv_hor,salt_adv_ver,         &
            avmols,salt_AH_method,salt_AH_const,salt_AH_Prt,  &
            salt_AH_stirr_const,sss_nudging_time,             &
            salt_check,min_salt,max_salt
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
!   integer, save :: Ncall = 0
!   Ncall = Ncall+1
!   write(debug,*) 'init_salinity() # ',Ncall
#endif

   LEVEL2 'init_salinity()'
   read(NAMLST,salt)

   if (avmols .lt. _ZERO_) then
      LEVEL3 'setting avmols to 0.'
      avmols = _ZERO_
   else
      LEVEL3 'avmols = ',real(avmols)
   end if

   if (.not. hotstart) then
      call init_salinity_field()
   end if

!  Sanity checks for advection specifications
   LEVEL3 'Advection of salinity'
   if (salt_adv_hor .eq. J7) stop 'init_salinity: J7 not implemented yet'
   call print_adv_settings_3d(salt_adv_split,salt_adv_hor,salt_adv_ver,_ZERO_)

   select case (salt_AH_method)
      case(0)
         LEVEL3 'salt_AH_method=0 -> horizontal diffusion of salt disabled'
      case(1)
         LEVEL3 'salt_AH_method=1 -> Using constant horizontal diffusivity of salt'
         if (salt_AH_const .lt. _ZERO_) then
              call getm_error("init_salinity()", &
                         "Constant horizontal diffusivity of salt <0");
         end if
         LEVEL4 real(salt_AH_const)
      case(2)
         LEVEL3 'salt_AH_method=2 -> using LES parameterisation'
         LEVEL4 'turbulent Prandtl number: ',real(salt_AH_Prt)
         deformC_3d =.true.
         deformX_3d =.true.
         deformUV_3d=.true.
         if (Am_method .eq. AM_LES) then
            les_mode = LES_BOTH
         else
            les_mode = LES_TRACER
         end if
      case(3)
         LEVEL3 'salt_AH_method=3 -> SGS stirring parameterisation'
         if (salt_AH_stirr_const .lt. _ZERO_) then
              call getm_error("init_salinity()", &
                         "salt_AH_stirr_const <0");
         end if
         LEVEL4 'stirring constant: ',real(salt_AH_stirr_const)
         deformC_3d =.true.
         deformX_3d =.true.
         deformUV_3d=.true.
         calc_stirr=.true.
      case default
         call getm_error("init_salinity()", &
                         "A non valid salt_AH_method has been chosen");
   end select

   if (metforcing .and. sss_nudging_time.gt._ZERO_) then
      nudge_sss = .True.
      LEVEL3 'nudging of SSS enabled'
      LEVEL4 'sss_nudging_time=',real(sss_nudging_time)
      allocate(sss(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_salinity: Error allocating memory (sss)'
      select case(met_method)
         case (METEO_CONST)
            if (sss_const .gt. _ZERO_) then
               LEVEL4 'constant sss=',real(sss_const)
               sss = sss_const
            else
               call getm_error("init_salinity()", &
                               "non-positive sss_const")
            end if
         case (METEO_FROMFILE,METEO_FROMEXT)
            LEVEL4 'sss read'
      end select
   end if

   LEVEL3 'salt_check=',salt_check
   if (salt_check .ne. 0) then
      LEVEL4 'doing sanity check on salinity'
      LEVEL4 'min_salt=',min_salt
      LEVEL4 'max_salt=',max_salt
      if (salt_check .gt. 0) then
         LEVEL4 'out-of-bound values result in termination of program'
      end if
      if (salt_check .lt. 0) then
         LEVEL4 'out-of-bound values result in warnings only'
      end if

      if (.not. hotstart) then

      call check_3d_fields(imin,jmin,imax,jmax,kmin,kmax,az, &
                           S,min_salt,max_salt,status)
      if (status .gt. 0) then
         if (salt_check .gt. 0) then
            call getm_error("init_salinity()", &
                            "out-of-bound values encountered")
         end if
         if (salt_check .lt. 0) then
            LEVEL1 'init_salinity(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if

      end if

   end if


#ifdef DEBUG
   write(debug,*) 'Leaving init_salinity()'
   write(debug,*)
#endif
   return
   end subroutine init_salinity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_salinity_field - initialisation of the salinity field
! \label{sec-init-salinity-field}
!
! !INTERFACE:
   subroutine init_salinity_field()
!
! !DESCRIPTION:
! Initialisation of the salinity field as specified by salt\_method
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
   integer                   :: i,j,k,n
   integer, parameter        :: nmax=10000
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: status
!EOP
!-------------------------------------------------------------------------
!BOC

   select case (salt_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value ',real(salt_const)
         forall(i=imin:imax,j=jmin:jmax, az(i,j) .ne. 0) &
                S(i,j,:) = salt_const
      case(2)
         LEVEL3 'using profile'
         LEVEL4 trim(salt_file)
         call read_profile(salt_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,kmax, &
                           az,H,hn,S)
      case(3)
         LEVEL3 'interpolating from 3D field'
         LEVEL4 trim(salt_file)
         call get_3d_field(salt_file,salt_name,salt_field_no,.true.,S)
      case default
         FATAL 'Not valid salt_method specified'
         stop 'init_salinity'
   end select

   S(:,:,0) = -9999.0
   forall(i=imin:imax,j=jmin:jmax, az(i,j).eq.0) S(i,j,:) = -9999.0

   call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,D_TAG)
   call wait_halo(D_TAG)
   call mirror_bdy_3d(S,D_TAG)

#ifdef DEBUG
   write(debug,*) 'Leaving init_salinity_field()'
   write(debug,*)
#endif
   return
   end subroutine init_salinity_field
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_salinity - salinity equation \label{sec-do-salinity}
!
! !INTERFACE:
   subroutine do_salinity(n)
!
! !DESCRIPTION:
!
! Here, one time step for the salinity equation is performed.
! First, preparations for the call to the advection schemes are
! made, i.e.\ calculating the necessary metric coefficients.
! After the call to the advection schemes, which actually perform
! the advection (and horizontal diffusion) step as an operational
! split step, the tri-diagonal matrix for calculating the
! new salinity by means of a semi-implicit central scheme for the
! vertical diffusion is set up.
! There are no source terms on the right hand sides.
! The subroutine is completed by solving the tri-diagonal linear
! equation by means of a tri-diagonal solver.
!
! Also here, there are some specific options for single test cases
! selected by compiler options.
!
! !USES:
   use advection_3d, only: do_advection_3d
   use variables_3d, only: dt,cnpar,hn,ho,nuh,uu,vv,ww,hun,hvn
   use domain,       only: imin,imax,jmin,jmax,kmax,az
   use getm_timers, only: tic,toc,TIM_SALT,TIM_SALTH,TIM_MIXANALYSIS
   use variables_3d, only: do_numerical_analyses_3d
   use variables_3d, only: nummix_S,nummix_S_old,nummix_S_int
   use variables_3d, only: phymix_S,phymix_S_int
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in) :: n
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,rc
   REALTYPE, POINTER         :: Res(:)
   REALTYPE, POINTER         :: auxn(:),auxo(:)
   REALTYPE, POINTER         :: a1(:),a2(:),a3(:),a4(:)
  REALTYPE                   :: S2(I3DFIELD)
  integer                    :: status
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_salinity() # ',Ncall
#endif
   call tic(TIM_SALT)

   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .eq. 1) then
            S(i,j,kmax) = S(i,j,kmax)*(_ONE_-fwf_int(i,j)/ho(i,j,kmax))
         end if
      end do
   end do

#ifdef _NUMERICAL_ANALYSES_OLD_
   if (do_numerical_analyses_3d) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
! OMP-note: The following array-based line could be implemented
!    with OMP as a WORKSHARE construct. However, it would require a dedicated
!    PARALLEL region (as the various advection schemes have their own regions),
!    so the overhead of the contruct would be rather large.
      S2 = S**2
      call do_advection_3d(dt,S2,uu,vv,ww,hun,hvn,ho,hn,                           &
                           salt_adv_split,salt_adv_hor,salt_adv_ver,_ZERO_,H_TAG)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if
#endif

   call do_advection_3d(dt,S,uu,vv,ww,hun,hvn,ho,hn,                           &
                        salt_adv_split,salt_adv_hor,salt_adv_ver,_ZERO_,H_TAG, &
                        nvd=nummix_S)

#ifdef _NUMERICAL_ANALYSES_OLD_
   if (do_numerical_analyses_3d) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
      call numerical_mixing(S2,S,nummix_S_old,nummix_S_int)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if
#endif

   if (salt_AH_method .gt. 0) then
!     S is not halo updated after advection
      call tic(TIM_SALTH)
      call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)
      call toc(TIM_SALTH)

      call tracer_diffusion(S,hn,salt_AH_method,salt_AH_const,salt_AH_Prt,salt_AH_stirr_const, &
                            phymix_S)
   end if

   if (do_numerical_analyses_3d) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
      call physical_mixing(S,avmols,phymix_S,phymix_S_int,salt_AH_method)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if


! OMP-NOTE: Pointers are used to for each thread to use its
!           own work storage.
!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,rc)                                             &
!$OMP    PRIVATE(Res,auxn,auxo,a1,a2,a3,a4)

! Each thread allocates its own HEAP storage:
   allocate(Res(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (Res)'
   allocate(auxn(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (auxn)'
   allocate(auxo(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (auxo)'
   allocate(a1(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (a1)'
   allocate(a2(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (a2)'
   allocate(a3(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (a3)'
   allocate(a4(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'do_salinity: Error allocating memory (auxo)'

! Note: We do not need to initialize these work arrays.
!   Tested BJB 2009-09-25.



!  Advection and vertical diffusion and of salinity
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            if (kmax.gt.1) then
!     Auxilury terms, old and new time level,
               do k=1,kmax-1
                  auxo(k)=_TWO_*(1-cnpar)*dt*(nuh(i,j,k)+avmols)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=_TWO_*   cnpar *dt*(nuh(i,j,k)+avmols)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
               end do

!        Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)
               a2(k)=hn(i,j,k)+auxn(k-1)
               a4(k)=S(i,j,k)*(hn(i,j,k)-auxo(k-1))+S(i,j,k-1)*auxo(k-1)
               if (nudge_sss) then
!                 implicit nudging
                  a2(kmax) = a2(kmax) + dry_z(i,j)*hn(i,j,kmax)*dt/sss_nudging_time
                  a4(kmax) = a4(kmax) + dry_z(i,j)*hn(i,j,kmax)*dt/sss_nudging_time*sss(i,j)
!                 explicit nudging
                  !a4(kmax) = a4(kmax) - dry_z(i,j)*dt*hn(i,j,kmax)*(S(i,j,kmax)-sss(i,j))/sss_nudging_time
               end if

!        Matrix elements for inner layers
               do k=2,kmax-1
                  a3(k)=-auxn(k  )
                  a1(k)=-auxn(k-1)
                  a2(k)=hn(i,j,k)+auxn(k)+auxn(k-1)
                  a4(k)=S(i,j,k+1)*auxo(k)                           &
                       +S(i,j,k  )*(hn(i,j,k)-auxo(k)-auxo(k-1))     &
                       +S(i,j,k-1)*auxo(k-1)
               end do

!        Matrix elements for bottom layer
               k=1
               a3(k)=-auxn(k  )
               a2(k)=hn(i,j,k)+auxn(k)
               a4(k)=S(i,j,k+1)*auxo(k)                              &
                    +S(i,j,k  )*(hn(i,j,k)-auxo(k))

               call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)

               do k=1,kmax
                  S(i,j,k)=Res(k)
               end do

            end if
         end if
      end do
   end do
!$OMP END DO

! Each thread must deallocate its own HEAP storage:
   deallocate(Res,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (Res)'
   deallocate(auxn,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (auxn)'
   deallocate(auxo,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (auxo)'
   deallocate(a1,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (a1)'
   deallocate(a2,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (a2)'
   deallocate(a3,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (a3)'
   deallocate(a4,stat=rc)
   if (rc /= 0) stop 'do_salinity: Error deallocating memory (a4)'

!$OMP END PARALLEL


   if (salt_check .ne. 0 .and. mod(n,abs(salt_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,kmin,kmax,az, &
                           S,min_salt,max_salt,status)
      if (status .gt. 0) then
         if (salt_check .gt. 0) then
            call getm_error("do_salinity()", &
                            "out-of-bound values encountered")
         end if
         if (salt_check .lt. 0) then
            LEVEL1 'do_salinity(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if

   call toc(TIM_SALT)
#ifdef DEBUG
   write(debug,*) 'Leaving do_salinity()'
   write(debug,*)
#endif
   return
   end subroutine do_salinity
!EOC

!-----------------------------------------------------------------------

   end module salinity

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
