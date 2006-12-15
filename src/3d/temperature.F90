!$Id: temperature.F90,v 1.20 2006-12-15 10:25:42 kbk Exp $
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
   use domain, only: imin,jmin,imax,jmax,H,az
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use variables_3d, only: T,hn,adv_schemes,kmin
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_temperature, do_temperature
!
! !PRIVATE DATA MEMBERS:
   integer                   :: temp_method=1,temp_format=2
   character(len=PATH_MAX)   :: temp_file="t_and_s.nc"
   character(len=32)         :: temp_name='temp'
   REALTYPE                  :: temp_const=20.
   integer                   :: temp_hor_adv=1,temp_ver_adv=1
   integer                   :: temp_adv_split=0
   REALTYPE                  :: temp_AH=-1.
   integer                   :: temp_check=0
   REALTYPE                  :: min_temp=-2.,max_temp=35.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
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
   subroutine init_temperature(adv_method)
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
   integer                   :: k,i,j,n
   integer, parameter        :: nmax=10000
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: temp_field_no=1
   NAMELIST /temp/ &
            temp_method,temp_const,temp_file,                 &
            temp_format,temp_name,temp_field_no,              &
            temp_hor_adv,temp_ver_adv,temp_adv_split,temp_AH, &
            temp_check,min_temp,max_temp
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_temperature() # ',Ncall
#endif

#ifdef NS_FRESHWATER_LENSE_TEST
temp_field_no=1
#endif
#ifdef MED_15X15MINS_TEST
temp_field_no=1
#endif

   LEVEL2 'init_temperature()'
   read(NAMLST,temp)

   select case (temp_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value'
         forall(i=iimin:iimax,j=jjmin:jjmax, az(i,j) .ne. 0) &
                T(i,j,:) = temp_const
      case(2)
         LEVEL3 'using profile'
         call read_profile(temp_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,az,H,       &
                           iimin,jjmin,iimax,jjmax,kmax,hn,T)
      case(3)
         LEVEL3 'interpolating from 3D field'
         call get_field(temp_file,temp_name,temp_field_no,T)
      case default
         FATAL 'Not valid temp_method specified'
         stop 'init_temperature'
   end select

!  Sanity checks for advection specifications
   LEVEL3 'temp_hor_adv=   ',temp_hor_adv
   LEVEL3 'temp_ver_adv=   ',temp_ver_adv
   LEVEL3 'temp_adv_split= ',temp_adv_split
   if(temp_hor_adv .eq. 1) then
      temp_adv_split=-1
      if(temp_ver_adv .ne. 1) then
         LEVEL3 "setting temp_ver_adv to 1 - since temp_hor_adv is 1"
         temp_ver_adv=1
      end if
   end if
   LEVEL3 "horizontal: ",trim(adv_schemes(temp_hor_adv))
   LEVEL3 "vertical:   ",trim(adv_schemes(temp_ver_adv))

   select case (temp_adv_split)
      case (-1)
      case (0)
         select case (temp_hor_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=0: temp_hor_adv not valid (2-6)")

         end select
         select case (temp_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=0: temp_ver_adv not valid (2-6)")
         end select
         LEVEL3 "1D split --> full u, full v, full w"
      case (1)
         select case (temp_hor_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=1: temp_hor_adv not valid (2-6)")
         end select
         select case (temp_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=1: temp_ver_adv not valid (2-6)")
         end select
         LEVEL3 "1D split --> half u, half v, full w, half v, half u"
      case (2)
         select case (temp_hor_adv)
            case (2,7)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=2: temp_hor_adv not valid (2,7)")
         end select
         select case (temp_ver_adv)
            case (2,3,4,5,6)
            case default
               call getm_error("init_3d()", &
                    "temp_adv_split=2: temp_ver_adv not valid (2-6)")
         end select
         LEVEL3 "2D-hor, 1D-vert split --> full uv, full w"
      case default
   end select

   call update_3d_halo(T,T,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)
   call wait_halo(D_TAG)
   call mirror_bdy_3d(T,D_TAG)

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
   use domain,       only: iimin,iimax,jjmin,jjmax,kmax,az,au,av
   use meteo,        only: swr,shf
   use parameters,   only: rho_0,cp
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use parameters, only: avmolt
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,rc
   REALTYPE                  :: Res(0:kmax)
   REALTYPE                  :: auxn(1:kmax-1),auxo(1:kmax-1)
   REALTYPE                  :: a1(0:kmax),a2(0:kmax)
   REALTYPE                  :: a3(0:kmax),a4(0:kmax)
   REALTYPE                  :: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE                  :: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE                  :: area_inv(I2DFIELD)
   REALTYPE                  :: swr_loc,shf_loc
   REALTYPE                  :: zz,rad(0:1000),A=0.58,g1=0.35,g2=23.0
   integer                   :: status
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_temperature() # ',Ncall
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
   call do_advection_3d(dt,T,uu,vv,ww,hun,hvn,ho,hn,   &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,     &
                        temp_hor_adv,temp_ver_adv,temp_adv_split,temp_AH)
#ifdef FRESHWATER_LENSE_TEST
   T(iimin:iimin+3,jjmin:jjmax,1:kmax)=10.
   T(iimax-3:iimax,jjmin:jjmax,1:kmax)=10.
   T(iimin:iimax,jjmin:jjmin+3,1:kmax)=10.
   T(iimin:iimax,jjmax-3:jjmax,1:kmax)=10.
#endif

!  Solar radiation and vertical diffusion of temperature

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1) then

! Solar radiation
            swr_loc=swr(i,j)
            shf_loc=shf(i,j)
            if (T(i,j,kmax).le.-0.0575*S(i,j,kmax)) then  ! use most primitive 
                                                          ! sea ice model ...
               shf_loc=max(_ZERO_,shf_loc)
            end if
            rad(kmax)=(swr_loc+shf_loc)/(rho_0*cp)
            zz = _ZERO_
            do k=kmax-1,0,-1
               zz=zz+hn(i,j,k+1)
               rad(k)=swr_loc/(rho_0*cp)*(A*exp(-zz/g1)+(1-A)*exp(-zz/g2))
            end do

            if (kmax.gt.1) then
!     Auxilury terms, old and new time level,
               do k=1,kmax-1
                  auxo(k)=2.*(1-cnpar)*dt*(nuh(i,j,k)+avmolt)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=2.*   cnpar *dt*(nuh(i,j,k)+avmolt)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
               end do

!        Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)
               a2(k)=hn(i,j,k)+auxn(k-1)
               a4(k)=T(i,j,k)*(hn(i,j,k)-auxo(k-1))+T(i,j,k-1)*auxo(k-1)  &
                     +dt*(rad(k)-rad(k-1))

!        Matrix elements for inner layers
               do k=2,kmax-1
                  a3(k)=-auxn(k  )
                  a1(k)=-auxn(k-1)
                  a2(k)=hn(i,j,k)+auxn(k)+auxn(k-1)
                  a4(k)=T(i,j,k+1)*auxo(k)                          &
                       +T(i,j,k  )*(hn(i,j,k)-auxo(k)-auxo(k-1))    &
                       +T(i,j,k-1)*auxo(k-1)                        &
                       +dt*(rad(k)-rad(k-1))
               end do

!        Matrix elements for bottom layer
               k=1
               a3(k)=-auxn(k  )
               a2(k)=hn(i,j,k)+auxn(k)
               a4(k)=T(i,j,k+1)*auxo(k)                           &
                    +T(i,j,k  )*(hn(i,j,k)-auxo(k))               &
                    +dt*(rad(k)-rad(k-1))

               call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)

               do k=1,kmax
                  T(i,j,k)=Res(k)
               end do

            end if
         end if
      end do
   end do

   if (temp_check .ne. 0 .and. mod(n,abs(temp_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,az,       &
                           iimin,jjmin,iimax,jjmax,kmax, &
                           kmin,T,min_temp,max_temp,status)
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
