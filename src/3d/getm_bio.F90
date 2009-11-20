!$Id: getm_bio.F90,v 1.4 2009-11-20 09:17:13 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: getm_bio()
!
! !INTERFACE:
   module getm_bio
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g,rho_0
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use time, only: secondsofday
   use variables_2d, only: D
   use variables_3d, only: taub
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d, only: nuh,T,S,rad,rho,light
   use variables_3d, only: cc3d
   use advection_3d, only: do_advection_3d
   use meteo, only: swr,u10,v10
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   use bio, only: init_bio, init_var_bio, set_env_bio, do_bio
   use bio, only: bio_calc
   use bio_var, only: numc
   use bio_var, only: cc,ws,var_names,var_units,var_long
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_getm_bio, do_getm_bio
   integer, public           :: bio_init_method=0
!
! !PRIVATE DATA MEMBERS:
   integer         :: bio_hor_adv=1
   integer         :: bio_ver_adv=1
   integer         :: bio_adv_split=1
   REALTYPE        :: bio_AH=-1.
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
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: getm_bio.F90,v $
!  Revision 1.4  2009-11-20 09:17:13  kb
!  updated GOTM/BIO interface
!
!  Revision 1.3  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.2  2007-05-14 12:16:22  kbk
!  added support for biology - via GOTM
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_bio
!
! !INTERFACE:
   subroutine init_getm_bio(nml_file)
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
   integer, parameter        :: unit_bio=63
   integer                   :: rc
   integer                   :: i,j,n
   character(len=PATH_MAX)   :: bio_init_file
   integer                   :: bio_init_format, bio_field_no

   namelist /getm_bio_nml/ bio_init_method, &
                           bio_init_file,bio_init_format,bio_field_no, &
                           bio_hor_adv,bio_ver_adv,bio_adv_split,bio_AH
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_bio()'

   call init_bio(NAMLST2,'bio.nml',unit_bio,kmax)

   if (bio_calc) then

      call init_var_bio

      allocate(cc3d(numc,I3DFIELD),stat=rc)         ! biological fields
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d)'
      cc3d = _ZERO_
!      cc3d(10,:,:,:) = 0.0001
      open(NAMLST2,status='unknown',file=trim(nml_file))
      read(NAMLST2,NML=getm_bio_nml)
      close(NAMLST2)

      LEVEL2 "settings related to 3D biological calculations"
      LEVEL3 'bio_hor_adv=   ',bio_hor_adv
      LEVEL3 'bio_ver_adv=   ',bio_ver_adv
      LEVEL3 'bio_adv_split= ',bio_adv_split
      LEVEL3 'bio_AH=        ',bio_AH

      select case (bio_init_method)
         case(0)
            LEVEL3 'getting initial bio fields from hotstart file'
         case(1)
            LEVEL3 "initial biological fields from namelist - bio_<model>.nml"
            do j=jmin,jmax
               do i=imin,imax
                  if (az(i,j) .ge. 1 ) then
                     cc3d(:,i,j,:)=cc
                  end if
               end do
            end do
         case(2)
            LEVEL3 'reading initial bio-fields from ',trim(bio_init_file)
            do n=1,numc
               LEVEL4 'inquiring ',trim(var_names(n))
               call get_field(bio_init_file,trim(var_names(n)),bio_field_no, &
                              cc3d(n,:,:,:))
            end do
         case default
            FATAL 'Not valid bio_init_method specified'
            stop 'init_getm_bio()'
      end select

      do n=1,numc
         call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az, &
                             imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      end do

#ifndef STATIC
      allocate(delxu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxu)'

      allocate(delxv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxv)'

      allocate(delyu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyu)'

      allocate(delyv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyv)'

      allocate(area_inv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (area_inv)'

      allocate(ff(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (ff)'
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
   end subroutine init_getm_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_getm_bio()
!
! !INTERFACE:
   subroutine do_getm_bio(dt)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer         :: n
   integer         :: i,j,k
   REALTYPE        :: h1d(0:kmax),T1d(0:kmax),S1d(0:kmax),rho1d(0:kmax)
   REALTYPE        :: nuh1d(0:kmax),rad1d(0:kmax),light1d(0:kmax)
   REALTYPE        :: bioshade1d(0:kmax)
   REALTYPE        :: wind_speed,I_0
!EOP
!-----------------------------------------------------------------------
!BOC

!  First we do all the vertical processes
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then
#ifdef GOTM_V3
            I_0=swr(i,j)
            h1d=hn(i,j,:)
            T1d=T(i,j,:)
            S1d=S(i,j,:)
            nuh1d=nuh(i,j,:)
            light1d=light(i,j,:)
            cc=cc3d(:,i,j,:)
            call do_bio(kmax,I_0,dt,h1d,T1d,S1d,nuh1d,light1d,bioshade1d)
            cc3d(:,i,j,:)=cc
            light(i,j,:)=bioshade1d
#else
            h1d=hn(i,j,:)
            T1d=T(i,j,:)
            S1d=S(i,j,:)
            rho1d=rho(i,j,:)
            nuh1d=nuh(i,j,:)
            rad1d=rad(i,j,:)
            if (allocated(u10) .and. allocated(v10)) then
               wind_speed=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
            else
               wind_speed=_ZERO_
            end if
            if (allocated(swr)) then
               I_0=swr(i,j)
            else
               I_0=_ZERO_
            end if
            light1d=light(i,j,:)
            cc=cc3d(:,i,j,:)
            call set_env_bio(kmax,dt,-D(i,j),sqrt(taub(i,j)), &
                             h1d,T1d,S1d,rho1d,nuh1d,rad1d,   &
                             wind_speed,I_0,secondsofday,w_adv_ctr_=0)
            call do_bio()
            cc3d(:,i,j,:)=cc
!            light(i,j,:)=bioshade1d
#endif
         end if
      end do
   end do

!  then we do the advection of the biological variables
   do n=1,numc

#if 1
      ff = cc3d(n,:,:,:)
      call do_advection_3d(dt,ff,uu,vv,ww,hun,hvn,ho,hn, &
              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
              bio_hor_adv,bio_ver_adv,bio_adv_split,bio_AH)

      call update_3d_halo(ff,ff,az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)
      cc3d(n,:,:,:) = ff
#else
      call do_advection_3d(dt,cc3d(n,:,:,:),uu,vv,ww,hun,hvn,ho,hn, &
              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
              bio_hor_adv,bio_ver_adv,bio_adv_split,bio_AH)

      call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
      call wait_halo(D_TAG)
#endif
   end do

   return
   end subroutine do_getm_bio
!EOC

!-----------------------------------------------------------------------

   end module getm_bio

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
