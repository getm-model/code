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
   use variables_2d, only: fwf_int
   use variables_3d, only: S,hn,adv_schemes,kmin
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
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
   REALTYPE                  :: salt_const=35.
   integer                   :: salt_hor_adv=1,salt_ver_adv=1
   integer                   :: salt_adv_split=0
   REALTYPE                  :: salt_AH=-_ONE_
   integer                   :: salt_check=0
   REALTYPE                  :: min_salt=0.,max_salt=40.
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
   subroutine init_salinity(adv_method)
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: adv_method
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
            salt_hor_adv,salt_ver_adv,salt_adv_split,salt_AH, &
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
   LEVEL3 'salt_hor_adv=   ',salt_hor_adv
   LEVEL3 'salt_ver_adv=   ',salt_ver_adv
   LEVEL3 'salt_adv_split= ',salt_adv_split
   if(salt_hor_adv .eq. 1) then
      salt_adv_split=-1
      if(salt_ver_adv .ne. 1) then
         LEVEL3 "setting salt_ver_adv to 1 - since salt_hor_adv is 1"
         salt_ver_adv=1
      end if
   end if
   LEVEL3 "horizontal: ",trim(adv_schemes(salt_hor_adv))
   LEVEL3 "vertical:   ",trim(adv_schemes(salt_ver_adv))

   select case (salt_adv_split)
      case (-1)
      case (0)
         select case (salt_hor_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=0: salt_hor_adv not valid (2-6)")
         end select
         select case (salt_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=0: salt_ver_adv not valid (2-6)")
         end select
         LEVEL3 "1D split --> full u, full v, full w"
      case (1)
         select case (salt_hor_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=1: salt_hor_adv not valid (2-6)")
         end select
         select case (salt_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=1: salt_ver_adv not valid (2-6)")
         end select
         LEVEL3 "1D split --> half u, half v, full w, half v, half u"
      case (2)
         select case (salt_hor_adv)
            case (2,7)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=2: salt_hor_adv not valid (2,7)")
         end select
         select case (salt_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "salt_adv_split=2: salt_ver_adv not valid (2-6)")
         end select
         LEVEL3 "2D-hor, 1D-vert split --> full uv, full w"
      case default
   end select

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
! \label{sec-init-salinity}
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
#ifdef PECS_TEST
   integer                   :: cc(1:30)
#endif
   integer, parameter        :: nmax=100
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: status
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef NS_FRESHWATER_LENSE_TEST
salt_field_no=1
#endif
#ifdef MED_15X15MINS_TEST
salt_field_no=1
#endif

   select case (salt_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value'
         forall(i=imin:imax,j=jmin:jmax, az(i,j) .ne. 0) &
                S(i,j,:) = salt_const
      case(2)
         LEVEL3 'using profile'
         call read_profile(salt_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,kmax, &
                           az,H,hn,S)
      case(3)
         LEVEL3 'interpolating from 3D field'
         call get_field(salt_file,salt_name,salt_field_no,S)
#ifdef SALTWEDGE_TEST
      case(4)
      !I need this here for hotstart of salinity!!!
         LEVEL3 'initializing with #ifdef SALTWEDGE'
         S =  _ZERO_
         do i=1,100
            do k=1,kmax
               S(i,2,k)=30.*(1.- tanh(float(i-1)*0.05))
            end do      
         end do   
#endif
      case default
         FATAL 'Not valid salt_method specified'
         stop 'init_salinity'
   end select

#ifdef ARKONA_TEST
   do i=100,135
      do j=256,257
         if (az(i,j).ge.1) S(i,j,0:kmax) = 25.
      end do
   end do
   do i=26,27
      do j=77,100
         S(i,j,0:kmax) = 8.
      end do
   end do
#endif
#ifdef INTERLEAVING_TEST
   S(2:6,2,1:20) = 12.
#endif
#ifdef SLOPE_TEST
   do i=81,82
      do j=42,43
      S(i,j,0:kmax) = 25.
      end do
   end do
#endif
#ifdef BALTIC_SLICE_TEST
  j=2
   if (imax.eq.102) then
      do i=2,21
        S(i,j,0:kmax) = 25.
      end do
   end if
   if (imax.eq.302) then
      do i=2,61
        S(i,j,0:kmax) = 25.
      end do
   end if
   if (imax.eq.902) then
      do i=2,181
        S(i,j,0:kmax) = 25.
      end do
   end if
#endif
!#else
!#ifdef PECS_TEST
!   S = 10.
!   do i=1,160
!      read(98,*) cc(1:30)
!      do j=1,30
!         if (cc(j).eq.1) then
!            S(i,j+1,1)=20.
!         end if
!      end do
!   end do
!#endif
!#endif

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
   use domain,       only: imin,imax,jmin,jmax,kmax,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use parameters, only: avmols
   use getm_timers, only: tic, toc, TIM_SALT, TIM_MIXANALYSIS
   use variables_3d, only: do_mixing_analysis
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
#ifdef SALTWEDGE_TEST
   REALTYPE                  :: SRelax,kk
#endif
  REALTYPE                   :: S2(I3DFIELD)
  REALTYPE                   :: delxu(I2DFIELD),delxv(I2DFIELD)
  REALTYPE                   :: delyu(I2DFIELD),delyv(I2DFIELD)
  REALTYPE                   :: area_inv(I2DFIELD)
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

   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            S(i,j,kmax) = S(i,j,kmax)*hn(i,j,kmax) &
                          /(hn(i,j,kmax)+fwf_int(i,j))
         end if
      end do
   end do

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

   if (do_mixing_analysis) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
! OMP-note: The following array-based line could be implemented
!    with OMP as a WORKSHARE construct. However, it would require a dedicated
!    PARALLEL region (as the various advection schemes have their own regions),
!    so the overhead of the contruct would be rather large.
      S2 = S**2
      call do_advection_3d(dt,S2,uu,vv,ww,hun,hvn,ho,hn,    &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,   &
                        salt_hor_adv,salt_ver_adv,salt_adv_split,salt_AH)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if

   call do_advection_3d(dt,S,uu,vv,ww,hun,hvn,ho,hn,    &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,   &
                        salt_hor_adv,salt_ver_adv,salt_adv_split,salt_AH)

   if (do_mixing_analysis) then
      call toc(TIM_SALT)
      call tic(TIM_MIXANALYSIS)
      call numerical_mixing(S2,S,nummix3d_S,nummix2d_S)
      call physical_mixing(S,avmols,phymix3d_S,phymix2d_S)
      call toc(TIM_MIXANALYSIS)
      call tic(TIM_SALT)
   end if

#ifdef PECS_TEST
   S(imin:imin,jmin:jmax,1:kmax)=10*_ONE_
   S(imax:imax,jmin:jmax,1:kmax)=10*_ONE_
#endif

#ifdef SALTWEDGE_TEST
   SRelax=30.
   j=2
   do k=1,kmax
      do i=1,100
         kk=  1.- tanh(float(i-1)*0.05)
         S(i,j,k)=(1.-kk)*S(i,j,k)+kk*SRelax
      end do
   end do
   S(imax-1,2,:)=0. !river
#endif


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
                  auxo(k)=2*(1-cnpar)*dt*(nuh(i,j,k)+avmols)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=2*   cnpar *dt*(nuh(i,j,k)+avmols)/ &
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

#ifdef ARKONA_TEST
   do i=100,135
      do j=256,257
         if (az(i,j).ge.1) S(i,j,0:kmax) = 25.
      end do
   end do
   do i=26,27
      do j=77,100
      S(i,j,0:kmax) = 8.
      end do
   end do
#endif

#ifdef SLOPE_TEST
   do i=81,82
      do j=42,43
      S(i,j,0:kmax) = 25.
      end do
   end do
#endif

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
