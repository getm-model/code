!$Id: temperature.F90,v 1.11 2005-04-25 07:55:50 kbk Exp $
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
!  Description still missing
!
! !USES:
   use exceptions
   use domain, only: imin,jmin,imax,jmax,H,az
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use variables_3d, only: T,hn,adv_schemes
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
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: temperature.F90,v $
!  Revision 1.11  2005-04-25 07:55:50  kbk
!  use more general frame for error handling - Umlauf
!
!  Revision 1.10  2004/01/06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!  Revision 1.9  2003/12/16 17:10:05  kbk
!  removed TABS
!
!  Revision 1.8  2003/12/16 16:13:51  kbk
!  forced ????_strang to 0 - needs clarification
!
!  Revision 1.7  2003/12/16 16:00:46  kbk
!  molecular diffusion for salt and temp (manuel)
!
!  Revision 1.6  2003/09/13 10:52:21  kbk
!  changed field_no to salt_field_no and temp_field_no
!
!  Revision 1.5  2003/08/03 08:13:09  kbk
!  added field_no to namelist
!
!  Revision 1.4  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:58  gotm
!  recovering after CVS crash
!
!  Revision 1.14  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.13  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.12  2001/08/31 15:35:48  bbh
!  nmax extended to 10000
!
!  Revision 1.11  2001/08/29 12:08:05  bbh
!  temp_adv and salt_adv needs to be public
!
!  Revision 1.10  2001/08/29 11:21:46  bbh
!  namelists read in salinity and temperature + initialisation
!
!  Revision 1.9  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.8  2001/07/26 13:14:06  bbh
!  Included adv_method in calls
!
!  Revision 1.7  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.6  2001/05/22 08:26:35  bbh
!  Added vertical diffusion
!
!  Revision 1.5  2001/05/18 12:01:57  bbh
!  Typo in #define LOCK_EXCHANGE_TEST
!
!  Revision 1.4  2001/05/18 09:51:07  bbh
!  Included az in call to update_3d_halo()
!
!  Revision 1.3  2001/05/18 09:45:16  bbh
!  Removed some LEVEL2 statements
!
!  Revision 1.2  2001/05/18 08:24:41  bbh
!  Advection of salinity and temperature
!
!  Revision 1.1  2001/05/03 20:20:33  bbh
!  Stubs for baroclinicity
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_temperature
!
! !INTERFACE:
   subroutine init_temperature(adv_method)
!
! !DESCRIPTION:
!  Description still missing
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
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: k,i,j,n
   integer, parameter        :: nmax=10000
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: temp_field_no=1
   NAMELIST /temp/ &
             temp_method,temp_const,temp_file,              &
             temp_format,temp_name,temp_field_no,           &
             temp_hor_adv,temp_ver_adv,temp_adv_split,temp_AH
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_temperature() # ',Ncall
#endif

#ifdef NS_NOMADS_TEST
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
   LEVEL3 "horizontal: ",trim(adv_schemes(temp_hor_adv))," of temperature"
   LEVEL3 "vertical:   ",trim(adv_schemes(temp_ver_adv))," of temperature"

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
         LEVEL2 "1D split --> full u, full v, full w"
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
         LEVEL2 "1D split --> half u, half v, full w, half v, half u"
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
         LEVEL2 "2D-hor, 1D-vert split --> full uv, full w"
      case default
   end select

   call update_3d_halo(T,T,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)
   call wait_halo(D_TAG)

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
! !IROUTINE:  do_temperature()
!
! !INTERFACE:
   subroutine do_temperature(n)
!
! !DESCRIPTION:
!  Description still missing
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
! !REVISION HISTORY:
!  See the log for the module
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
#ifdef NOMADS_TEST
   T(iimin:iimin,jjmin:jjmax,1:kmax)=10.
   T(iimax:iimax,jjmin:jjmax,1:kmax)=10.
   T(iimin:iimax,jjmin:jjmin,1:kmax)=10.
   T(iimin:iimax,jjmax:jjmax,1:kmax)=10.
#endif

!  Solar radiation and vertical diffusion of temperature

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1) then

! Solar radiation
            swr_loc=swr(i,j)
            shf_loc=shf(i,j)
            if (T(i,j,kmax).le.-0.0575*S(i,j,kmax)) then  ! use most primitive ice model ...
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

   call update_3d_halo(T,T,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)
   call wait_halo(D_TAG)

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
