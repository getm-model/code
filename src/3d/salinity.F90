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
   use domain, only: imin,jmin,imax,jmax,kmax,ioff,joff
   use domain, only: H,az
!KB   use get_field, only: get_3d_field
   use variables_2d, only: fwf_int
   use variables_3d, only: rk,S,hn,kmin
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG,H_TAG
   IMPLICIT NONE
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
   REALTYPE                  :: salt_AH=-_ONE_
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
   subroutine init_salinity()
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
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   integer                   :: status
   NAMELIST /salt/                                            &
            salt_method,salt_const,salt_file,                 &
            salt_format,salt_name,salt_field_no,              &
            salt_adv_split,salt_adv_hor,salt_adv_ver,salt_AH, &
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

   call init_salinity_field()

!  Sanity checks for advection specifications
   LEVEL3 'Advection of salinity'
   if (salt_adv_hor .eq. J7) stop 'init_salinity: J7 not implemented yet'
   call print_adv_settings_3d(salt_adv_split,salt_adv_hor,salt_adv_ver,salt_AH)

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

      if (salt_method .ne. 0) then
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

   S(:,:,0) = -9999._rk
   forall(i=imin:imax,j=jmin:jmax, az(i,j).eq.0) S(i,j,:) = -9999._rk

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
   use parameters, only: avmols
   use getm_timers, only: tic, toc, TIM_SALT, TIM_MIXANALYSIS
   use variables_3d, only: do_numerical_analyses
   use variables_3d, only: nummix3d_S,nummix2d_S
   use variables_3d, only: phymix3d_S,phymix2d_S
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
! Developers note:
!  The parentheses are set to minimize truncation errors for fwf_int=0
            S(i,j,kmax) = S(i,j,kmax)*            &
                          ( hn(i,j,kmax) / (hn(i,j,kmax)+fwf_int(i,j)) )
         end if
      end do
   end do

   if (do_numerical_analyses) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
! OMP-note: The following array-based line could be implemented
!    with OMP as a WORKSHARE construct. However, it would require a dedicated
!    PARALLEL region (as the various advection schemes have their own regions),
!    so the overhead of the contruct would be rather large.
      S2 = S**2
      call do_advection_3d(dt,S2,uu,vv,ww,hun,hvn,ho,hn,                           &
                           salt_adv_split,salt_adv_hor,salt_adv_ver,salt_AH,H_TAG)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if

   call do_advection_3d(dt,S,uu,vv,ww,hun,hvn,ho,hn,                            &
                        salt_adv_split,salt_adv_hor,salt_adv_ver,salt_AH,H_TAG)

   if (do_numerical_analyses) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)

      call numerical_mixing(S2,S,nummix3d_S,nummix2d_S)

      call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)
      call physical_mixing(S,salt_AH,avmols,phymix3d_S,phymix2d_S)

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
