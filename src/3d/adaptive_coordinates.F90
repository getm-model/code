!$Id: adaptive_coordinates.F90,v 1.6 2010-02-23 08:23:34 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  adaptive vertical coordinates
! \label{sec-adaptive-coordinates}
!
! !INTERFACE:
   subroutine adaptive_coordinates(first,hotstart)
!
! !DESCRIPTION:
! The vertical grid adaptivity is partially given by a vertical diffusion
! equation for the vertical layer positions, with diffusivities being
! proportional to shear, stratification and distance from the boundaries.
! In the horizontal, the grid can be smoothed with respect to $z$-levels,
! grid layer slope and density. Lagrangian tendency of the grid movement
! is supported. The adaptive terrain-following grid can be set to be an
! Eulerian-Lagrangian grid, a hybrid $\sigma$-$\rho$ or $\sigma$-$z$ grid
! and combinations of these with great flexibility. With this, internal
! flow structures such as thermoclines can be well resolved and
! followed by the grid. A set of idealised examples is presented in
! Hofmeister et al. (2009), which show that the introduced adaptive grid
! strategy reduces pressure gradient errors and numerical mixing significantly.
!
! For the configuration of parameters, a seperate namelist file adaptcoord.inp
! has to be given with parameters as following:
! \\
! faclag - Factor on Lagrangian coords., 0.le.faclag.le.1\\
! facdif - Factor on thickness filter,   0.le.faclag.le.1\\
! fachor - Factor on position filter,  0.le.faclag.le.1\\
! cNN - dependence on stratification\\
! cSS - dependence on shear\\
! cdd - dep. on distance from surface and bottom\\
! Snorm - Typical value of absolute shear\\
! NNnorm - Typical value of BVF squared\\
! dsurf - reference value for surface/bottom distance [m]\\
! tgrid - Time scale of grid adaptation [s]\\
! \\
! The parameters cNN,cSS,cdd are used for the vertical adaption and 
! have to be less or equal 1 in sum. The difference to 1 is describing a
! background value which forces the coordinates back to a sigma distribution.
! The values ddu and ddl from the domain namelist are used for weighting
! the zooming to surface and bottom if cdd>0.
! The compiler option PREADAPT=number allows for a pre-adaption of
! coordinates to the initial density field and bathymetry. The number
! defines the number of macro-timesteps used for the preadaption.
! The initial temperature and salinity fields are re-interpolated
! onto the adapted grid afterwards.
!
! !USES:
   use domain, only: ga,imin,imax,jmin,jmax,kmax,H,HU,HV,az,au,av
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hvo,hun,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   use variables_3d, only: kmin_pmz,kumin_pmz,kvmin_pmz

! ADAPTIVE-BEGIN
   use  parameters,  only: g,rho_0
   use variables_3d, only: uu,vv,NN,SS
   use variables_3d, only: rho
   use domain,       only: ddu,ddl
   use halo_zones, only: update_3d_halo,wait_halo
   use halo_zones, only: H_TAG,U_TAG,V_TAG
#if defined CURVILINEAR || defined SPHERICAL
   use domain,       only: dxv,dyu,arcd1
#else
   use domain,       only: dx,dy,ard1
#endif

!ADAPTIVE-END

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
   logical, intent(in)                 :: hotstart
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister and Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc
   REALTYPE        :: kmaxm1
   REALTYPE        :: deltaiso
   REALTYPE, save, dimension(:),     allocatable  :: be
   REALTYPE, save, dimension(:),     allocatable  :: NNloc ! local NN vector
   REALTYPE, save, dimension(:),     allocatable  :: SSloc ! local SS vector
   REALTYPE, save, dimension(:),     allocatable  :: gaa   ! new relative coord.
   REALTYPE, save, dimension(:),     allocatable  :: gaaold! old relative coord.
   REALTYPE, save, dimension(:),     allocatable  :: aav ! total grid diffus.
   REALTYPE, save, dimension(:),     allocatable  :: avn ! NN-rel. grid diffus.
   REALTYPE, save, dimension(:),     allocatable  :: avs ! SS-rel. grid diffus.
   REALTYPE, save, dimension(:),     allocatable  :: avd ! dist.-rel. grid diff.
   REALTYPE, save, dimension(:,:,:), allocatable  :: zpos ! new pos. of z-levels
   REALTYPE, save, dimension(:,:,:), allocatable  :: zposo! old pos. of z-levels
   REALTYPE, save, dimension(:,:,:), allocatable  :: work2,work3
   REALTYPE     :: faclag=0.5  ! Factor on Lagrangian coords., 0.le.faclag.le.1
   REALTYPE     :: facdif=0.3  ! Factor on thickness filter,   0.le.faclag.le.1
   REALTYPE     :: fachor=0.7  ! Factor on position filter,  0.le.faclag.le.1
   REALTYPE     :: faciso=0.e-6 ! Factor for isopycnal tendency
   REALTYPE     :: depthmin=0.1
   REALTYPE     :: Ncrit=1.e-6
   integer      :: mhor=1 ! this number is experimental - it has to be 1 for now-
   integer      :: iw=2 ! stencil for isopycnal tendency
   REALTYPE     :: rm
   INTEGER      :: im,iii,jjj,ii
   integer      :: split=1     !splits the vertical adaption into #split steps
   REALTYPE     :: c1ad=0.4     ! dependence on NN
   REALTYPE     :: c2ad=0.3     ! dependence on SS
   REALTYPE     :: c3ad=0.3    ! distance from surface and bottom
   REALTYPE     :: c4ad=0.0    ! background value
   REALTYPE     :: Snorm=0.01  ! Typical value of absolute shear 
   REALTYPE     :: NNnorm=0.05  ! Typical value of BVF squared
   REALTYPE     :: dsurf=1.     ! reference value for surface/bottom distance
   REALTYPE     :: tgrid=86400. ! Time scale of grid adaptation
   REALTYPE     :: dtgrid
   REALTYPE     :: aau(0:kmax),bu(0:kmax)
   REALTYPE     :: cu(0:kmax),du(0:kmax)
   REALTYPE     :: facupper=1.0
   REALTYPE     :: faclower=1.0
   REALTYPE     :: cNN,cSS,cdd,csum
   REALTYPE     :: cbg=1.0

   integer,save :: count=0     
      namelist /adapt_coord/   faclag,facdif,fachor,faciso,  &
                        depthmin,Ncrit, &
                        cNN,cSS,cdd,cbg,Snorm,NNnorm, &
                        dsurf,tgrid,split
#if (defined PARALLEL && defined INPUT_DIR)
   character(len=PATH_MAX)   :: input_dir=INPUT_DIR
#else
   character(len=PATH_MAX)   :: input_dir='./'
#endif

   
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

#ifdef DEBUG
STDERR 'adaptive_coordinates()'
#endif

   if (first) then

!read namelist
      open(ADAPTNML,status='unknown',file=(trim(input_dir) // 'adaptcoord.inp'))
      read(ADAPTNML,adapt_coord)
      close(ADAPTNML)

      if (.not. allocated(ga)) allocate(ga(0:kmax),stat=rc)
         if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      do k=0,kmax
         ga(k) = k
      end do

      allocate(zpos(I3DFIELD),stat=rc)  ! z-coord. of interface
      if (rc /= 0) stop 'coordinates: Error allocating memory (zpos)'
      allocate(zposo(I3DFIELD),stat=rc)  ! z-coord. of interface
      if (rc /= 0) stop 'coordinates: Error allocating memory (zposo)'
      allocate(work2(I3DFIELD),stat=rc)  !
      if (rc /= 0) stop 'coordinates: Error allocating memory (work2)'
      allocate(work3(I3DFIELD),stat=rc)  !
      if (rc /= 0) stop 'coordinates: Error allocating memory (work3)'
      allocate(be(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (be)'
      allocate(NNloc(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (NNloc)'
      allocate(SSloc(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (SSloc)'
      allocate(avn(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (avn)'
      allocate(avs(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (avs)'
      allocate(avd(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (avd)'
      allocate(aav(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (aav)'
      allocate(gaa(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (gaa)'
      allocate(gaaold(0:kmax),stat=rc)     ! working space
      if (rc /= 0) STOP 'coordinates: Error allocating (gaaold)'
      kmaxm1= _ONE_/float(kmax)

      if (.not.hotstart) then
! Dirty way to read initial distribution (as equidistant sigma coordinates):
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
               hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
            end do
         end do
         do j=jmin-HALO,jmax+HALO
            do i=imin-1-HALO,imax+HALO
               huo(i,j,:)=(ssuo(i,j)+HU(i,j))*kmaxm1
               hun(i,j,:)=(ssun(i,j)+HU(i,j))*kmaxm1
            end do
         end do
         do j=jmin-1-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               hvo(i,j,:)=(ssvo(i,j)+HV(i,j))*kmaxm1
               hvn(i,j,:)=(ssvn(i,j)+HV(i,j))*kmaxm1
            end do
         end do

      else !hotstart
         ho=hn
      end if

      kmin=1
      kumin=1
      kvmin=1
      kmin_pmz=1
      kumin_pmz=1
      kvmin_pmz=1

!   factors for diffusion at surface and bottom
      faclower=max(_ZERO_,ddl)/(max(_ZERO_,ddu)+max(_ZERO_,ddl))
      facupper=max(_ZERO_,ddu)/(max(_ZERO_,ddu)+max(_ZERO_,ddl))
!   norm factors for diffusive adaption
      c1ad=max(_ZERO_,cNN)
      c2ad=max(_ZERO_,cSS)
      c3ad=max(_ZERO_,cdd)
      if (cbg .lt. 1.0) then
         c4ad=max(_ZERO_,cbg)
      else
         c4ad=max(_ZERO_,(1.0-c1ad-c2ad-c3ad))
      end if
      csum=c1ad+c2ad+c3ad+c4ad
         if (csum .gt. _ZERO_) then
            c1ad=c1ad/csum
            c2ad=c2ad/csum
            c3ad=c3ad/csum
         end if

   else !if not first

      ho=hn
! Lagrangian and thickness filtering step 
      do k=1,kmax
         do j=jmin+1-HALO,jmax-1+HALO
            do i=imin+1-HALO,imax-1+HALO
               hn(i,j,k)=ho(i,j,k)                       &
                  -((uu(i,j,k)*DYU-uu(i-1,j  ,k)*DYUIM1) &
                  +(vv(i,j,k)*DXV-vv(i  ,j-1,k)*DXVJM1)) &
                  *ARCD1*dt*faclag                       &
                  +(                                     &
                   (ho(i+1,j  ,k)-ho(i,j,k))*min(1,au(i  ,j  )) &
                  +(ho(i-1,j  ,k)-ho(i,j,k))*min(1,au(i-1,j  )) &
                  +(ho(i  ,j+1,k)-ho(i,j,k))*min(1,av(i  ,j  )) &
                  + (ho(i  ,j-1,k)-ho(i,j,k))*min(1,av(i ,j-1)) &
                  )*facdif*0.25

! Ensure smooth transition to cut-off around depth
               hn(i,j,k)=max(hn(i,j,k),depthmin)
            end do
         end do
      end do
      call hcheck(hn,ssen,H)

!  Update the halo zones
      call update_3d_halo(hn,hn,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

! First step provided Lagrangian modified hn with additional horizontal
! diffusion of hn. Consistency with total depth is ensured

! Horizontal diffusion of zpos can be repeated mhor times
! mhor should be 1, because NN has to be recalculated
! NN is used without interpolation onto the new grid, the local linear density
! profile assumes NN=const.

      do ii=1,mhor
! Prepare zpos field
         call htoz(hn,zpos)

! For problem zones (boundaries, thin layers), put zero 'horizontal diffusion'
         do k=1,kmax-1
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  if (au(i,j).ge.1) work2(i,j,k)=_ONE_
               end do
            end do
         end do
         do k=2,kmax
            do j=jmin-HALO,jmax+HALO
               do i=imin+1-HALO,imax+HALO
                  if ((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                     work2(i  ,j,k  )=0
                     work2(i  ,j,k-1)=0
                     work2(i-1,j,k  )=0
                     work2(i-1,j,k-1)=0
                  endif
               end do
            end do
         end do
         do k=1,kmax
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  if (av(i,j).ge.1) work3(i,j,k)=_ONE_
               end do
            end do
         end do
         do k=2,kmax
            do j=jmin+1-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                     work3(i,j,k)=0
                     work3(i,j,k-1)=0
                     work3(i,j-1,k)=0
                     work3(i,j-1,k-1)=0
                  endif
               end do
            end do
         end do

         zposo=zpos
         do k=1,kmax-1
            do j=jmin-HALO+2,jmax+HALO-2
               do i=imin-HALO+2,imax+HALO-2
                  if (faciso .ge. 0.00001) then
                     rm=0
                     im=0
                     do iii=max(imin-HALO,i-iw),min(imax+HALO,i+iw)
                        do jjj=max(jmin-HALO,j-iw),min(jmax+HALO,j+iw)
                           rm=rm+min(1,az(iii,jjj))* &
                              (rho(iii,jjj,k+1)+rho(iii,jjj,k))
                           im=im+min(1,az(iii,jjj))
                        end do
                     end do
                     deltaiso= 0.5* (rm/im-rho(i,j,k+1)-rho(i,j,k))/     &
                        ((-1.0)*rho_0/9.81*max(NN(i,j,k),Ncrit))*faciso
                     if (deltaiso.ge.hn(i,j,k+1)) &
                        deltaiso=hn(i,j,k+1)-depthmin
                     if (deltaiso.le.(-hn(i,j,k))) &
                        deltaiso=-hn(i,j,k)+depthmin 
                     else
                        deltaiso=0.0
                     end if
  
                     zpos(i,j,k)=zposo(i,j,k)+ &
                        (   &
                        (zposo(i+1,j,k)-zposo(i,j,k))*work2(i,j,k) - &
                        (zposo(i,j,k)-zposo(i-1,j,k))*work2(i-1,j,k)  &
                        + (zposo(i,j+1,k)-zposo(i,j,k))*work3(i,j,k) - &
                        (zposo(i,j,k)-zposo(i,j-1,k))*work3(i,j-1,k)  &
                        )*0.25*fachor    &
                        + deltaiso
               end do
            end do
         end do
! Local consistency
         call ztoh(zpos,hn,depthmin)
         do k=1,kmax
            do j=jmin-HALO+1,jmax+HALO-1
               do i=imin-HALO+1,imax+HALO-1
                  if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                     zpos(i,j,k)=zpos(i,j,k-1)+depthmin
                  endif
               end do
            end do
         end do
 
         call ztoh(zpos,hn,depthmin)
         call hcheck(hn,ssen,H)


      end do ! End of Horizontal diffusion of zpos repeated mhor times
! After horizontal diffusion updated new and depth consistent hn is available

      call htoz(hn,zpos)
! get old z positions (same as NN, SS)
      call htoz(ho,zposo)

      dtgrid=dt/float(split)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .ge. 1) then
               NNloc=NN(i,j,:) 
               SSloc=SS(i,j,:)
               do k=0,kmax
                  gaa(k)=(zpos(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
                  gaaold(k)=(zposo(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
               end do 
               do ii=1,split

                  call col_interpol(kmax-1,1,gaaold,NN(i,j,:),kmax-1,gaa,NNloc)
                  call col_interpol(kmax-1,1,gaaold,SS(i,j,:),kmax-1,gaa,SSloc)

!     Stratification
                  NNloc(kmax)=NNloc(kmax-1)       
                  NNloc(0)=NNloc(1)       
                  SSloc(kmax)=SSloc(kmax-1)       
                  SSloc(0)=SSloc(1)       
                  do k=1,kmax
                     avn(k)=min(_ONE_,max(_ZERO_,0.5*(NNloc(k)+NNloc(k-1)))/g &
                           *rho_0/NNnorm)
                  end do

!     Shear
                  do k=1,kmax
                     avs(k)=min(_ONE_,sqrt(max(_ZERO_,0.5*  &
                            (SSloc(k)+SSloc(k-1))))/Snorm)
                  end do

!     Distance from surface and bottom
                  do k=1,kmax
                     avd(k)=facupper*1./(dsurf-0.5*(gaa(k-1)+gaa(k))*H(i,j))+         &
                           faclower*1./(dsurf+(1.+0.5*(gaa(k-1)+gaa(k)))*H(i,j))
                  end do

!     Calculation of grid diffusivity
                  do k=1,kmax
                     aav(k)=H(i,j)/tgrid*(c1ad*avn(k)+c2ad*avs(k)+ &
                           c3ad*avd(k)+c4ad/H(i,j))
                     aav(k)=aav(k)*dtgrid*kmax**2/100.
!        Minimum layer thickness
                     if ((gaa(k)-gaa(k-1)).lt.0.001/float(kmax)) aav(k)=0.
                  end do

                  do k=1,kmax-1
                     aau(k)=-aav(k)
                     cu(k)=-aav(k+1)
                     bu(k)=1.-aau(k)-cu(k)
                     du(k)=gaa(k)
                  end do
                  cu(0)=0
                  bu(0)=1.
                  du(0)=-1.

                  bu(kmax)=1.
                  aau(kmax)=0.
                  du(kmax)=0.

                  call getm_tridiagonal(kmax,0,kmax,aau,bu,cu,du,gaa)
                 
               end do !split
               zpos(i,j,:)=gaa*(ssen(i,j)+H(i,j))+ssen(i,j)
            end if
         end do
      end do

! Back to hn
      call ztoh(zpos,hn,depthmin)
! Normally if positive defined vertical diffusion no check required

   end if !first

   call hcheck(hn,ssen,H)
! Finally derive interface grid sizes for uu and vv
! Interface treatment and check

!  Update the halo zones
   call update_3d_halo(hn,hn,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)
   call mirror_bdy_3d(hn,H_TAG)

! uu
   huo=hun
   do k=1,kmax
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            hun(i,j,k)=0.5*(hn(i,j,k)+hn(i+1,j,k))
         end do
      end do
   end do
   call hcheck(hun,ssun,HU)

! vv
   hvo=hvn
   do k=1,kmax
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            hvn(i,j,k)=0.5*(hn(i,j,k)+hn(i,j+1,k))
         end do
      end do
   end do
   call hcheck(hvn,ssvn,HV)

!  Update the halo zones for hun
   call update_3d_halo(hun,hun,au,imin,jmin,imax,jmax,kmax,U_TAG)
   call wait_halo(U_TAG)
   call mirror_bdy_3d(hun,U_TAG)
!  Update the halo zones for hvn
   call update_3d_halo(hvn,hvn,av,imin,jmin,imax,jmax,kmax,V_TAG)
   call wait_halo(V_TAG)
   call mirror_bdy_3d(hvn,V_TAG)

   if (first) then
      huo=hun
      hvo=hvn
      ho=hn
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving adaptive_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine adaptive_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
