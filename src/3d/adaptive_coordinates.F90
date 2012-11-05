#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adaptive vertical coordinates
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
! d\_vel - Typical velocity difference for scaling cNN adaption\\
! d\_dens - Typical density difference for scaling cSS adaption\\
! dsurf - reference value for surface/bottom distance [m]\\
! tgrid - Time scale of grid adaptation [s]\\
! preadapt - number of iterations for pre-adaptation\\
! \\
! The parameters cNN,cSS,cdd are used for the vertical adaption and
! have to be less or equal 1 in sum. The difference to 1 is describing a
! background value which forces the coordinates back to a sigma distribution.
! The values ddu and ddl from the domain namelist are used for weighting
! the zooming to surface and bottom if cdd>0.
! The option preadapt allows for a pre-adaption of
! coordinates to the initial density field and bathymetry. The number
! defines the number of iterations (change coordinates, vertically advect
! tracer, calculate vertical gradients)  used for the preadaption.
! The initial temperature and salinity fields are re-interpolated
! onto the adapted grid afterwards.
!
! !USES:
   use domain, only: ga,imin,imax,jmin,jmax,kmax,H,HU,HV,az,au,av
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hvo,hun,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   use variables_3d, only: kmin_pmz,kumin_pmz,kvmin_pmz
   use variables_3d, only: preadapt
   use m3d, only: calc_salt, calc_temp, bdy3d
   use vertical_coordinates,only: restart_with_ho,restart_with_hn

! ADAPTIVE-BEGIN
   use  parameters,  only: g,rho_0
   use variables_3d, only: uu,vv,SS
#ifndef NO_BAROCLINIC
   use variables_3d, only: NN
   use variables_3d, only: rho
#endif
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
! !OUTPUT PARAMETERS:
!   integer, intent(out)                :: preadapt
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
   REALTYPE     :: faclag=_ZERO_    ! Factor on Lagrangian coords., 0.le.faclag.le.1
   REALTYPE     :: facdif=3*_TENTH_ ! Factor on thickness filter,   0.le.faclag.le.1
   REALTYPE     :: fachor=_TENTH_   ! Factor on position filter,  0.le.faclag.le.1
   REALTYPE     :: faciso=_ZERO_    ! Factor for isopycnal tendency
   REALTYPE     :: depthmin=_ONE_/5
   REALTYPE     :: Ncrit=_ONE_/1000000
   integer      :: mhor=1    ! this number is experimental - it has to be 1 for now-
   integer      :: iw=2      ! stencil for isopycnal tendency
   REALTYPE     :: rm
   INTEGER      :: im,iii,jjj,ii
   integer      :: split=1           ! splits the vertical adaption into #split steps
   REALTYPE     :: c1ad=_ONE_/5      ! dependence on NN
   REALTYPE     :: c2ad=_ZERO_       ! dependence on SS
   REALTYPE     :: c3ad=_ONE_/5      ! distance from surface and bottom
   REALTYPE     :: c4ad=6*_TENTH_    ! background value
   REALTYPE     :: d_vel=_TENTH_     ! Typical value of absolute shear
   REALTYPE     :: d_dens=_HALF_     ! Typical value of BVF squared
   REALTYPE     :: dsurf=20*_ONE_    ! reference value for surface/bottom distance
   REALTYPE     :: tgrid=21600*_ONE_ ! Time scale of grid adaptation
   REALTYPE     :: dtgrid
   REALTYPE     :: aau(0:kmax),bu(0:kmax)
   REALTYPE     :: cu(0:kmax),du(0:kmax)
   REALTYPE     :: facupper=_ONE_
   REALTYPE     :: faclower=_ONE_
   REALTYPE     :: cNN,cSS,cdd,csum
   REALTYPE     :: cbg=6*_TENTH_
   REALTYPE     :: tfac_hor=3600*_ONE_ ! factor introducing a hor. adaption timescale = dt/tgrid_hor where tgrid_hor=tgrid/N
   integer      :: iip

   integer,save :: count=0
      namelist /adapt_coord/   faclag,facdif,fachor,faciso,  &
                        depthmin,Ncrit, &
                        cNN,cSS,cdd,cbg,d_vel,d_dens, &
                        dsurf,tgrid,split,preadapt
#if (defined GETM_PARALLEL && defined INPUT_DIR)
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

      if (.not. allocated(ga)) then
         allocate(ga(0:kmax),stat=rc)
         if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      end if
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
      LEVEL2 "allocated memory"

! In the case of no-baroclinic, the nnloc will never be updated,
! so we initiialize it here. Same for avn: Both just remain zero.
      NNloc(:) = _ZERO_
      avn(:)   = _ZERO_
      if (.not. restart_with_hn) then
         if (hotstart) then
            LEVEL2 'WARNING: assume equidistant sigma coordinates for hn'
         end if
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
            end do
         end do
      end if
!     only for backward compatibility
      if (hotstart) then
         if (.not. restart_with_ho) then
            LEVEL2 'WARNING: assume ho=hn'
            ho=hn
         end if
      else
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
            end do
         end do
      end if

      kmin=1
      kumin=1
      kvmin=1
      kmin_pmz=1
      kumin_pmz=1
      kvmin_pmz=1

!     factors for diffusion at surface and bottom
      faclower=max(_ZERO_,ddl)/(max(_ZERO_,ddu)+max(_ZERO_,ddl))
      facupper=max(_ZERO_,ddu)/(max(_ZERO_,ddu)+max(_ZERO_,ddl))
!   norm factors for diffusive adaption
#ifdef NO_BAROCLINIC
      c1ad=_ZERO_
#else
      c1ad=max(_ZERO_,cNN)
#endif
      c2ad=max(_ZERO_,cSS)
      c3ad=max(_ZERO_,cdd)
      if (cbg .lt. _ONE_) then
         c4ad=max(_ZERO_,cbg)
      else
         c4ad=max(_ZERO_,(_ONE_-c1ad-c2ad-c3ad))
      end if
      csum=c1ad+c2ad+c3ad+c4ad
      if (csum .gt. _ONE_) then
         c1ad=c1ad/csum
         c2ad=c2ad/csum
         c3ad=c3ad/csum
      end if

      ! possibly useful:
      !tfac_hor=dt*kmax/tgrid
      !LEVEL2 "horizontal time scale:",tfac_hor
      !
      ! used in Hofmeister et al. 2010:
      tfac_hor=_ONE_

   else !not first
      ho=hn
!     Lagrangian and thickness filtering step
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
                  +(ho(i  ,j-1,k)-ho(i,j,k))*min(1,av(i  ,j-1)) &
                  )*facdif*_QUART_ *tfac_hor

!              Ensure smooth transition to cut-off around depth
               hn(i,j,k)=max(hn(i,j,k),depthmin)
            end do
         end do
      end do
      call hcheck(hn,ssen,H)

!     Update the halo zones
      call update_3d_halo(hn,hn,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

! First step provided Lagrangian modified hn with additional horizontal
! diffusion of hn. Consistency with total depth is ensured

! Horizontal diffusion of zpos can be repeated mhor times
! mhor should be 1, because NN has to be recalculated
! NN is used without interpolation onto the new grid, the local linear density
! profile assumes NN=const.

      do ii=1,mhor
!        Prepare zpos field
         call htoz(hn,zpos)

!         For problem zones (boundaries, thin layers), put zero 'horizontal diffusion'
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
                     work2(i  ,j,k  )=_ZERO_
                     work2(i  ,j,k-1)=_ZERO_
                     work2(i-1,j,k  )=_ZERO_
                     work2(i-1,j,k-1)=_ZERO_
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
                     work3(i,j  ,k)   = _ZERO_
                     work3(i,j  ,k-1) = _ZERO_
                     work3(i,j-1,k)   = _ZERO_
                     work3(i,j-1,k-1) = _ZERO_
                  endif
               end do
            end do
         end do

         zposo=zpos
         do k=1,kmax-1
            do j=jmin-HALO+2,jmax+HALO-2
               do i=imin-HALO+2,imax+HALO-2
                  deltaiso=_ZERO_
#ifndef NO_BAROCLINIC
                  if (faciso .ge. _ONE_/100000) then
                     rm=_ZERO_
                     im=0
                     do iii=max(imin-HALO,i-iw),min(imax+HALO,i+iw)
                        do jjj=max(jmin-HALO,j-iw),min(jmax+HALO,j+iw)
                           rm=rm+min(1,az(iii,jjj))* &
                              (rho(iii,jjj,k+1)+rho(iii,jjj,k))
                           im=im+min(1,az(iii,jjj))
                        end do
                     end do
                     deltaiso= _HALF_* (rm/im-rho(i,j,k+1)-rho(i,j,k))/     &
                        ((-_ONE_)*rho_0/g*max(NN(i,j,k),Ncrit))* &
                        faciso * tfac_hor
                     if (deltaiso.ge.hn(i,j,k+1)) &
                        deltaiso=hn(i,j,k+1)-depthmin
                     if (deltaiso.le.(-hn(i,j,k))) &
                        deltaiso=-hn(i,j,k)+depthmin
                  end if
#endif
                  zpos(i,j,k)= zposo(i,j,k) + deltaiso                       &
                        +_QUART_*fachor *tfac_hor*(_ZERO_                    &
                          +(zposo(i+1,j  ,k)-zposo(i  ,j,k))*work2(i  ,j,k)  &
                          -(zposo(i,  j  ,k)-zposo(i-1,j,k))*work2(i-1,j,k)  &
                          +(zposo(i,  j+1,k)-zposo(i,  j,k))*work3(i,  j,k)  &
                          -(zposo(i,  j,  k)-zposo(i,j-1,k))*work3(i,j-1,k) )
               end do
            end do
         end do
!        Local consistency
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

!     After horizontal diffusion updated new and depth consistent hn is available
      call htoz(hn,zpos)
!     Get old z positions (same as NN, SS)
      call htoz(ho,zposo)

      dtgrid=dt/float(split)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .ge. 1) then
#ifndef NO_BAROCLINIC
               NNloc=NN(i,j,:)
#endif
               SSloc=SS(i,j,:)
               do k=0,kmax
                  gaa(k)   =( zpos(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
                  gaaold(k)=(zposo(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
               end do
               do ii=1,split
#ifndef NO_BAROCLINIC
!                 Stratification
                  call col_interpol(kmax-1,1,gaaold,NN(i,j,:),kmax-1,gaa,NNloc)
                  NNloc(kmax)=NNloc(kmax-1)
                  NNloc(0)=NNloc(1)
                  do k=1,kmax
                     avn(k)=min(_ONE_,max(_ZERO_,_HALF_*(NNloc(k)+NNloc(k-1)))/g &
                           *rho_0/d_dens)
                  end do
#endif
!                 Shear
                  call col_interpol(kmax-1,1,gaaold,SS(i,j,:),kmax-1,gaa,SSloc)
                  SSloc(kmax)=SSloc(kmax-1)
                  SSloc(0)=SSloc(1)
                  do k=1,kmax
                     avs(k)=min(_ONE_,sqrt(max(_ZERO_,_HALF_*  &
                            (SSloc(k)+SSloc(k-1))))/d_vel)
                  end do

!                 Distance from surface and bottom
                  do k=1,kmax
                     avd(k)=facupper*_ONE_/(dsurf-_HALF_*(gaa(k-1)+gaa(k))*H(i,j))+         &
                           faclower*_ONE_/(dsurf+(_ONE_+_HALF_*(gaa(k-1)+gaa(k)))*H(i,j))
                  end do

!                 Calculation of grid diffusivity
                  do k=1,kmax
                     aav(k)=H(i,j)/tgrid*(c1ad*avn(k)+c2ad*avs(k)+ &
                           c3ad*avd(k)+c4ad/H(i,j))
                     aav(k)=aav(k)*dtgrid*kmax**2/100.
!                    Minimum layer thickness
                     if ((gaa(k)-gaa(k-1)).lt.0.001/float(kmax)) aav(k)=0.
!                      if (hn(i,j,k).le.depthmin) aav(k)=0.
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

!     Back to hn
      call ztoh(zpos,hn,depthmin)
!     Normally if positive defined vertical diffusion no check required

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
            hun(i,j,k)=_QUART_*(ho(i,j,k)+ho(i+1,j,k)+hn(i,j,k)+hn(i+1,j,k))
         end do
      end do
   end do
!  KK-TODO: although the layer heights in the center points are consistent
!           with the total water depth, in the present implementation we
!           cannot rely on depth-coinciding layer heights in velocity points
   call hcheck(hun,ssun,HU)

! vv
   hvo=hvn
   do k=1,kmax
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            hvn(i,j,k)=_QUART_*(ho(i,j,k)+ho(i,j+1,k)+hn(i,j,k)+hn(i,j+1,k))
         end do
      end do
   end do
!  KK-TODO: although the layer heights in the center points are consistent
!           with the total water depth, in the present implementation we
!           cannot rely on depth-coinciding layer heights in velocity points
   call hcheck(hvn,ssvn,HV)

!  Update the halo zones for hun
   call update_3d_halo(hun,hun,au,imin,jmin,imax,jmax,kmax,U_TAG)
   call wait_halo(U_TAG)
   call mirror_bdy_3d(hun,U_TAG)
!  Update the halo zones for hvn
   call update_3d_halo(hvn,hvn,av,imin,jmin,imax,jmax,kmax,V_TAG)
   call wait_halo(V_TAG)
   call mirror_bdy_3d(hvn,V_TAG)

!  only for backward compatibility
   if (first) then
      huo=hun
      hvo=hvn
      if (.not. hotstart) then
         ho=hn
      end if
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
