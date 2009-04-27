!$Id: adaptive_coordinates.F90,v 1.5 2009-04-27 07:36:08 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  adaptive vertical coordinates
! \label{sec-adaptive-coordinates}
!
! !INTERFACE:
   subroutine adaptive_coordinates(first)
!
! !DESCRIPTION:
!  For Richard to do
!
! !USES:
#if 0
   use domain, only: ga,imin,imax,jmin,jmax,kmax,H,HU,HV,az,au,av
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hvo,hun,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   use variables_3d, only: kmin_pmz,kumin_pmz,kvmin_pmz
! ADAPTIVE-BEGIN
   use  parameters,  only: g,rho_0
   use variables_3d, only: uu,vv,NN,SS
   use variables_3d, only: rho
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
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister and Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc
   REALTYPE        :: kmaxm1
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
   REALTYPE     :: faclag=0.0  ! Factor on Lagrangian coords., 0.le.faclag.le.1
   REALTYPE     :: facdif=0.0  ! Factor on thickness filter,   0.le.faclag.le.1
   REALTYPE     :: fachor=1.0  ! Factor on position filter,  0.le.faclag.le.1
   REALTYPE     :: faciso=0.e-6
   REALTYPE     :: depthmin=0.1
   REALTYPE     :: Ncrit=1.e-6
   integer      :: mhor=1
   integer      :: iw=1 ! iw=15
   REALTYPE     :: rm
   INTEGER      :: im,iii,jjj,ii
   integer      :: split=5
   REALTYPE     :: c1ad=0.3     ! dependence on NN
   REALTYPE     :: c2ad=0.2     ! dependence on SS
   REALTYPE     :: c3ad=0.3     ! distance from surface and bottom
   REALTYPE     :: c4ad=0.2     ! background value
   REALTYPE     :: Snorm=0.01  ! Typical value of absolute shear 
   REALTYPE     :: NNnorm=0.05  ! Typical value of BVF squared
   REALTYPE     :: dsurf=1.     ! reference value for surface/bottom distance
   REALTYPE     :: tgrid=86400. ! Time scale of grid adaptation
   REALTYPE     :: dtgrid
   REALTYPE     :: aau(0:kmax),bu(0:kmax)
   REALTYPE     :: cu(0:kmax),du(0:kmax)
   integer,save :: count=0     
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

STDERR 'adaptive_coordinates()'

   if (first) then
      if (.not. allocated(ga)) allocate(ga(0:kmax),stat=rc)
         if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      do k=0,kmax
         ga(k) = k
      end do

            write(*,*) 'begin first '
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

! Dirty way to read initial distribution (as equidistant sigma coordinates):
            do j=jmin,jmax
               do i=imin,imax
                  ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
                  hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
               end do
            end do
            do j=jmin,jmax
               do i=imin-1,imax
                  huo(i,j,:)=(ssuo(i,j)+HU(i,j))*kmaxm1
                  hun(i,j,:)=(ssun(i,j)+HU(i,j))*kmaxm1
               end do
            end do
            do j=jmin-1,jmax
               do i=imin,imax
                  hvo(i,j,:)=(ssvo(i,j)+HV(i,j))*kmaxm1
                  hvn(i,j,:)=(ssvn(i,j)+HV(i,j))*kmaxm1
               end do
            end do
            kmin=1
            kumin=1
            kvmin=1
            kmin_pmz=1
            kumin_pmz=1
            kvmin_pmz=1
   end if ! first


! DX and DY scaling should be included....
         ho=hn
! Lagrangian and thickness filtering step 
         do k=1,kmax
            do j=jmin+1,jmax-1
               do i=imin+1,imax-1
                  hn(i,j,k)=ho(i,j,k)                       &
                     -((uu(i,j,k)*DYU-uu(i-1,j  ,k)*DYUIM1) &
                     +(vv(i,j,k)*DXV-vv(i  ,j-1,k)*DXVJM1)) &
                     *ARCD1*dt*faclag                       &
                     +(                                     &
                      (ho(i+1,j  ,k)-ho(i,j,k))*au(i  ,j  ) &
                     +(ho(i-1,j  ,k)-ho(i,j,k))*au(i-1,j  ) &
                     +(ho(i  ,j+1,k)-ho(i,j,k))*av(i  ,j  ) &
                     + (ho(i  ,j-1,k)-ho(i,j,k))*av(i ,j-1) &
                     )*facdif*0.25

! Ensure smooth transition to cut-off around depth
                  hn(i,j,k)=max(hn(i,j,k),depthmin)
               end do
            end do
         end do

         call hcheck(hn,ssen,h)
! First step provided Lagrangian modified hn with additional horizontal
! diffusion of hn. Consistency with total depth is ensured

! Horizontal diffusion of zpos repeated mhor times
         do ii=1,mhor
! Prepare zpos field
         call htoz(hn,zpos)

! For problem zones (boundaries, thin layers), put zero 'horizontal diffusion'
         do k=1,kmax-1
            do j=jmin,jmax
               do i=imin,imax
                  if (au(i,j).eq.1) work2(i,j,k)=_ONE_
               end do
            end do
         end do
         do k=2,kmax
            do j=jmin,jmax
               do i=imin+1,imax
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
            do j=jmin,jmax
               do i=imin,imax
                  if (av(i,j).eq.1) work3(i,j,k)=_ONE_
               end do
            end do
         end do
         do k=2,kmax
            do j=jmin+1,jmax
               do i=imin,imax
                  if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                     work3(i,j,k)=0
                     work3(i,j,k-1)=0
                     work3(i,j-1,k)=0
                     work3(i,j-1,k-1)=0
                  endif
               end do
            end do
         end do
! Check for  bounds on i j etc and possible dirichlet conditions on z

! Dirty BC
           zposo=zpos
              do k=1,kmax-1
                 do j=jmin+1,jmax-1
                    do i=imin+1,imax-1
                       rm=0
                       im=0
                       do iii=max(imin,i-iw),min(imax,i+iw)
                          do jjj=max(jmin,j-iw),min(jmax,j+iw)
                             rm=rm+az(iii,jjj)*(rho(iii,jjj,k+1)+rho(iii,jjj,k))
                             im=im+az(iii,jjj)
                          end do
                       end do

                       zpos(i,j,k)=zposo(i,j,k)+ &
                          (   &
                          (zposo(i+1,j,k)-zposo(i,j,k))*work2(i,j,k) - &
                          (zposo(i,j,k)-zposo(i-1,j,k))*work2(i-1,j,k)  &
                          + (zposo(i,j+1,k)-zposo(i,j,k))*work3(i,j,k) - &
                          (zposo(i,j,k)-zposo(i,j-1,k))*work3(i,j-1,k)  &
                           )*0.25*fachor    &
                           +( &
                           0.5* (rm/im-rho(i,j,k+1)-rho(i,j,k))/     &
                           (max(NN(i,j,k),Ncrit)))*faciso
                    end do
                 end do
              end do
! Local consistency
              do k=1,kmax
                 do j=jmin,jmax
                    do i=imin,imax
                       if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                          zpos(i,j,k)=zpos(i,j,k-1)+depthmin
                       endif
                    end do
                 end do
              end do
 
              call ztoh(zpos,hn,depthmin)
              call hcheck(hn,ssen,h)

           end do ! End of Horizontal diffusion of zpos repeated mhor times
! After horizontal diffusion updated new and depth consistent hn is available

        call htoz(hn,zpos)
        dtgrid=dt/float(split)
        do j=jmin,jmax
           do i=imin,imax
              if (az(i,j) .eq. 1) then
              NNloc=NN(i,j,:) 
              SSloc=SS(i,j,:) 
              do k=0,kmax
                 gaa(k)=(zpos(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
                 gaaold(k)=gaa(k)
              end do 
              do ii=1,split
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
                    avd(k)=1./(dsurf-0.5*(gaa(k-1)+gaa(k))*H(i,j))+         &
                           1./(dsurf+(1.+0.5*(gaa(k-1)+gaa(k)))*H(i,j))
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

                 call col_interpol(kmax-1,1,gaaold,NN(i,j,:),kmax-1,gaa,NNloc)
                 call col_interpol(kmax-1,1,gaaold,SS(i,j,:),kmax-1,gaa,SSloc)
       
              end do 
              zpos(i,j,:)=gaa*(ssen(i,j)+H(i,j))+ssen(i,j)
              end if
           end do
        end do


! Back to hn
        call ztoh(zpos,hn,depthmin)
! Normally if positive defined vertical diffusion no check required

        hn=ho+0.1*(hn-ho) 
        call hcheck(hn,ssen,h)

! Finally derive interface grid sizes for uu and vv
! Interface treatment and check
! uu
        huo=hun
        do k=1,kmax
           do j=jmin,jmax
              do i=imin,imax-1
                 hun(i,j,k)=0.5*(hn(i,j,k)+hn(i+1,j,k))
              end do
           end do
        end do
! maybe not allowed in imax....
        call hcheck(hun,ssun,hu)
! vv
        hvo=hvn
           do k=1,kmax
              do j=jmin,jmax-1
                 do i=imin,imax
                    hvn(i,j,k)=0.5*(hn(i,j,k)+hn(i,j+1,k))
                 end do
              end do
           end do
! maybe not allowed in jmax....
           call hcheck(hvn,ssvn,hv)

#ifdef DEBUG
   write(debug,*) 'Leaving adaptive_coordinates()'
   write(debug,*)
#endif

#endif
   return
   end subroutine adaptive_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
