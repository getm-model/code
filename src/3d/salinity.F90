!$Id: salinity.F90,v 1.14 2005-04-25 07:55:50 kbk Exp $
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
!  Description still missing
!
! !USES:
   use exceptions
   use domain, only: imin,jmin,imax,jmax,ioff,joff
#ifdef HAIDVOGEL_TEST
   use domain, only: iextr,jextr
#endif
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use domain, only: H,az
   use variables_3d, only: S,hn,adv_schemes
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_salinity, do_salinity
!
! !PRIVATE DATA MEMBERS:
   integer                   :: salt_method=1,salt_format=2
   character(len=PATH_MAX)   :: salt_file="t_and_s.nc"
   character(len=32)         :: salt_name='salt'
   REALTYPE                  :: salt_const=35.
   integer                   :: salt_hor_adv=1,salt_ver_adv=1
   integer                   :: salt_adv_split=0
   REALTYPE                  :: salt_AH=-_ONE_
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: salinity.F90,v $
!  Revision 1.14  2005-04-25 07:55:50  kbk
!  use more general frame for error handling - Umlauf
!
!  Revision 1.13  2004/08/06 15:16:12  hb
!  Compiler option ARKONA_TEST included
!
!  Revision 1.12  2004/07/29 19:46:32  hb
!  For compiler option NOMADS_TEST: some lines shortened
!
!  Revision 1.11  2004/04/21 09:14:03  lars
!  added ifdef for baltic_slice
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
!  Revision 1.12  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.11  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.10  2001/08/29 12:08:05  bbh
!  temp_adv and salt_adv needs to be public
!
!  Revision 1.9  2001/08/29 11:21:46  bbh
!  namelists read in salinity and temperature + initialisation
!
!  Revision 1.8  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.7  2001/07/26 13:14:06  bbh
!  Included adv_method in calls
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/22 08:26:12  bbh
!  Fixed vertical diffusion
!
!  Revision 1.4  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.3  2001/05/18 09:51:07  bbh
!  Included az in call to update_3d_halo()
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
! !IROUTINE: init_salinity
!
! !INTERFACE:
   subroutine init_salinity(adv_method)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
#ifdef NOMADS_TEST
   use domain, only: dx
#endif
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
  integer                    :: i,j,k,n
#ifdef PECS_TEST
   integer                   :: cc(1:30)
#endif
#ifdef NOMADS_TEST
   INTEGER                   :: ic,jc
   REALTYPE                  :: dist
#endif
   integer, parameter        :: nmax=100
   REALTYPE                  :: zlev(nmax),prof(nmax)
   integer                   :: salt_field_no=1
   NAMELIST /salt/                                          &
            salt_method,salt_const,salt_file,               &
            salt_format,salt_name,salt_field_no,            &
            salt_hor_adv,salt_ver_adv,salt_adv_split,salt_AH
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_salinity() # ',Ncall
#endif

#ifdef NS_NOMADS_TEST
salt_field_no=1
#endif
#ifdef MED_15X15MINS_TEST
salt_field_no=1
#endif

   LEVEL2 'init_salinity()'
   read(NAMLST,salt)

   select case (salt_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value'
         forall(i=iimin:iimax,j=jjmin:jjmax, az(i,j) .ne. 0) &
                S(i,j,:) = salt_const
      case(2)
         LEVEL3 'using profile'
         call read_profile(salt_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,az,H,  &
                           iimin,jjmin,iimax,jjmax,kmax,hn,S)
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
   LEVEL3 "horizontal: ",trim(adv_schemes(salt_hor_adv))," of salinity"
   LEVEL3 "vertical:   ",trim(adv_schemes(salt_ver_adv))," of salinity"

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

#ifdef NOMADS_TEST
   S=34.85
   ic=nint(iimax/2.)
   jc=nint(jjmax/2.)
   do i=1,iimax
      do j=1,jjmax
         do k=kmax/2,kmax
            dist=sqrt((float(i)-float(ic))**2+(float(j)-float(jc))**2)*dx
            if (dist.le.3000.) then
               S(i,j,k)=1.1*(dist/(3000.))**8+33.75
            else
               S(i,j,k)=34.85
            end if
         end do
      end do
   end do
#endif
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
#ifdef BALTIC_SLICE_TEST
  j=2
   if (iimax.eq.102) then
      do i=2,21
        S(i,j,0:kmax) = 25.
      end do
   end if
   if (iimax.eq.302) then
      do i=2,61
        S(i,j,0:kmax) = 25.
      end do
   end if
   if (iimax.eq.902) then
      do i=2,181
        S(i,j,0:kmax) = 25.
      end do
   end if
#endif
#ifdef HAIDVOGEL_TEST
   STDERR 'HAIDVOGEL: salinity= ',iimin,iimax,i+ioff,iextr/2
   do i=iimin-1,iimax+1
      if(i+ioff .le. iextr/2) then
         S(i,jjmin-1:jjmax+1,0:kmax) = 6.4102564
         S(i,jjmin-1:jjmax+1,0:kmax) = 5.
      else
         S(i,jjmin-1:jjmax+1,0:kmax) = 0.
      end if
   end do
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

   call update_3d_halo(S,S,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)
   call wait_halo(D_TAG)

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
! !IROUTINE:  do_salinity()
!
! !INTERFACE:
   subroutine do_salinity(n)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   use advection_3d, only: do_advection_3d
   use variables_3d, only: dt,cnpar,hn,ho,nuh,uu,vv,ww,hun,hvn
   use domain,       only: iimin,iimax,jjmin,jjmax,kmax,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use parameters, only: avmols
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
#ifdef NOMADS_TEST
   REALTYPE                  :: SRelax,kk
#endif
#ifdef SALTWEDGE_TEST
   REALTYPE                  :: SRelax,kk
#endif
  REALTYPE                   :: delxu(I2DFIELD),delxv(I2DFIELD)
  REALTYPE                   :: delyu(I2DFIELD),delyv(I2DFIELD)
  REALTYPE                   :: area_inv(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_salinity() # ',Ncall
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
   call do_advection_3d(dt,S,uu,vv,ww,hun,hvn,ho,hn,    &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,   &
                        salt_hor_adv,salt_ver_adv,salt_adv_split,salt_AH)

#ifdef PECS_TEST
   S(iimin:iimin,jjmin:jjmax,1:kmax)=10.
   S(iimax:iimax,jjmin:jjmax,1:kmax)=10.
#endif

#ifdef NOMADS_TEST
   SRelax=34.85
   do i=1,iimax
      do j=1,jjmax
         do k=1,kmax
            if (az(i,j).eq.2) S(i,j,k)=SRelax
            kk=0.5625
            if ((((i.eq.2).or.(i.eq.iimax-1)).and.(j.ge.2).and.   &
                (j.le.jjmax-1)).or.(((j.eq.2).or.(j.eq.jjmax-1))  &
                .and.(i.ge.2).and.(i.le.iimax-1)))                &
                S(i,j,k)=(1.-kk)*S(i,j,k)+kk*SRelax
            kk=0.25
            if ((((i.eq.3).or.(i.eq.iimax-2)).and.(j.ge.3).and.   &
                (j.le.jjmax-2)).or.(((j.eq.3).or.(j.eq.jjmax-2))  &
                .and.(i.ge.3).and.(i.le.iimax-2)))                &
                S(i,j,k)=(1.-kk)*S(i,j,k)+kk*SRelax
            kk=0.0625
            if ((((i.eq.4).or.(i.eq.iimax-3)).and.(j.ge.4).and.   &
                (j.le.jjmax-3)).or.(((j.eq.4).or.(j.eq.jjmax-3))  &
                .and.(i.ge.4).and.(i.le.iimax-3)))                &
                S(i,j,k)=(1.-kk)*S(i,j,k)+kk*SRelax
         end do
      end do
  end do
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

!  Advection and vertical diffusion and of salinity

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1) then
            if (kmax.gt.1) then
!     Auxilury terms, old and new time level,
               do k=1,kmax-1
                  auxo(k)=2.*(1-cnpar)*dt*(nuh(i,j,k)+avmols)/ &
                             (hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=2.*   cnpar *dt*(nuh(i,j,k)+avmols)/ &
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

   call update_3d_halo(S,S,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)
   call wait_halo(D_TAG)

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
