!$Id: coordinates.F90,v 1.6 2004-04-21 15:18:15 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  coordinates() - defines the vertical coordinate system.
!
! !INTERFACE:
   subroutine coordinates(cord_type,cord_relax,maxdepth)
!
! !DESCRIPTION:
!
! !USES:
   use  parameters, only: g
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG,HU_TAG,HV_TAG
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,H,HU,HV,az,au,av,min_depth
   use domain, only: ga,ddu,ddl,d_gamma,gamma_surf
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hun,hvo,hvn
   use variables_3d, only: kmin_pmz,kumin_pmz,kvmin_pmz
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
!JMB
   use variables_3d, only: uu,vv,NN,SS
   use variables_3d, only: rho
   use domain, only: dx,dy,ard1
! to test ...
   use variables_3d, only: S

!JMB END

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: cord_type
   REALTYPE, intent(in)                :: cord_relax
   REALTYPE, intent(in)                :: maxdepth
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: coordinates.F90,v $
!  Revision 1.6  2004-04-21 15:18:15  hb
!  Algorithms for adaptive grids added
!
!  Revision 1.5  2004/01/05 13:23:27  kbk
!  Poor Man's Z-coordinates
!
!  Revision 1.4  2003/09/02 14:45:46  kbk
!  calculate in HALO-zones instead of using update_3d_halo()
!
!  Revision 1.3  2003/04/23 12:16:27  kbk
!  added calls to wait_halo()
!
!  Revision 1.2  2003/04/07 16:27:32  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:53  gotm
!  recovering after CVS crash
!
!  Revision 1.15  2001/10/26 07:42:27  bbh
!  Correct values for sigma and general cordinates in ga
!
!  Revision 1.14  2001/10/23 14:15:55  bbh
!  Moved ga from coordinates.F90 to domain.F90
!
!  Revision 1.13  2001/10/23 12:43:48  bbh
!  Forgot to calculate initial values hn, hun, hvn - when general vertical coodinates
!
!  Revision 1.12  2001/10/22 09:26:41  bbh
!  Added cord_relax
!
!  Revision 1.11  2001/10/17 07:43:59  bbh
!  Relaxation of layer depth + cleaning
!
!  Revision 1.10  2001/09/14 09:01:58  bbh
!  Re-ordered
!
!  Revision 1.9  2001/09/01 17:10:25  bbh
!  Vertical coordinate definition now specified via namelist
!
!  Revision 1.8  2001/08/31 15:44:44  bbh
!  general vertical coordinates added
!
!  Revision 1.7  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.6  2001/05/18 10:00:50  bbh
!  Added masks to calls to update_3d_halo()
!
!  Revision 1.5  2001/05/18 08:15:49  bbh
!  Cosmetics
!
!  Revision 1.4  2001/05/16 07:02:49  bbh
!  Finished equidistant and non-equidistant sigma coordinates
!
!  Revision 1.3  2001/05/15 11:45:38  bbh
!  Added zooming
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc,kk
   REALTYPE        :: tmp,kmaxm1,alpha
   REALTYPE        :: HH,zz,r,hnmin=0.01
   logical, save   :: first=.true.,equiv_sigma=.false.
   logical         :: kminset
   REALTYPE, save, dimension(:),     allocatable  :: dga,be,sig,zlev
   REALTYPE, save, dimension(:,:,:), allocatable  :: gga,zzz
!JMB
   REALTYPE, save, dimension(:),     allocatable  :: NNloc,SSloc
   REALTYPE, save, dimension(:),     allocatable  :: gaa,gaaold
   REALTYPE, save, dimension(:),     allocatable  :: aav,av1,av2,av3
   REALTYPE, save, dimension(:),     allocatable  :: jmb1,jmb2,jmb3,jmb4
   REALTYPE, save, dimension(:,:,:), allocatable  :: zpos
   REALTYPE, save, dimension(:,:,:), allocatable  :: zposo
   REALTYPE, save, dimension(:,:,:), allocatable  :: FF
   REALTYPE, save, dimension(:,:,:), allocatable  :: work,work2,work3
   REALTYPE     :: dtm1
   REALTYPE     :: FACLAG=0.
   REALTYPE     :: FACDIF=0.1
   REALTYPE     :: FACVER=0.095
   REALTYPE     :: FACHOR=0.00
   REALTYPE     :: FACISO=0.5
   REALTYPE     :: depthmin=0.05
   REALTYPE     :: SSnorm=2.0
   REALTYPE     :: NNnorm=0.2
   integer      :: MHOR=0
   integer      :: MVER=15
   integer      :: IW=15
   REALTYPE     :: RM,RLIM1,RLIM2,RLIM3,RLIM4
   INTEGER      :: IM,III,JJJ,II
   integer      :: split=5
   REALTYPE     :: c1ad=0.4,c2ad=0.5,c3ad=0.1,c4ad=0.0,dsurf=2.,Tgrid=17280.
   REALTYPE     :: dtgrid
   REALTYPE     :: aau(0:kmax),bu(0:kmax)
   REALTYPE     :: cu(0:kmax),du(0:kmax)
!JMB END
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

   if (first) then
      first = .false.
      if (.not. allocated(ga)) allocate(ga(0:kmax),stat=rc)
      if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      select case (cord_type)
         case (1) ! sigma coordinates
            if (ddu .le. _ZERO_ .and. ddl .le. _ZERO_) then
               equiv_sigma=.true.
               ga(0) = -_ONE_
               do k=1,kmax
                  ga(k) = ga(k-1) + _ONE_/kmax
               end do
               ga(kmax) = _ZERO_
            else
               ! Non-equidistant sigma coordinates
               ! This zooming routine is from Antoine Garapon, ICCH, DK
               if (ddu .lt. _ZERO_) ddu=_ZERO_
               if (ddl .lt. _ZERO_) ddl=_ZERO_
               allocate(dga(0:kmax),stat=rc)
               if (rc /= 0) STOP 'coordinates: Error allocating (dga)'
               ga(0)= -_ONE_
               dga(0)= _ZERO_
               do k=1,kmax
                  ga(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
                  ga(k)=ga(k)/(tanh(ddl)+tanh(ddu)) - _ONE_
                  dga(k)=ga(k)-ga(k-1)
               end do
            end if
            kmin=1
            kumin=1
            kvmin=1
            kmin_pmz=1
            kumin_pmz=1
            kvmin_pmz=1
         case (2) ! z-level
            allocate(zlev(0:kmax),stat=rc)
            if (rc /= 0) stop 'coordinates: Error allocating (zlev)'
            allocate(gga(I3DFIELD),stat=rc)  ! dimensionless gamma-coordinate
            if (rc /= 0) stop 'coordinates: Error allocating memory (gga)'
            ! Calculate profile of zlevels
            zlev(0)=-maxdepth
            do k=1,kmax
               zlev(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
               zlev(k)=zlev(k)/(tanh(ddl)+tanh(ddu)) - _ONE_
               zlev(k)=zlev(k)*maxdepth
            end do
            kmin_pmz=1
            do j=jjmin,jjmax
               do i=iimin,iimax
                 HH=max(sseo(i,j)+H(i,j),min_depth)
                 kminset=.true.
                 if (az(i,j).eq.0) then
                    do k=1,kmax
                       hn(i,j,k)=HH/float(kmax)
                       gga(i,j,k)=hn(i,j,k)/HH
                    end do
                    kmin_pmz(i,j)=1
                  else
                     k=kmax
111                     if (zlev(k-1)-(k-1)*hnmin.ge.-H(i,j)) then
                          hn(i,j,k)=zlev(k)-zlev(k-1)
                    else
                       hn(i,j,k)=max(zlev(k)-(-H(i,j)+(k-1)*hnmin),hnmin)
                       if (kminset.and.hn(i,j,k).eq.hnmin) then
                          kmin_pmz(i,j)=k+1
                          kminset=.false.
                       end if
                    end if
                    gga(i,j,k)=hn(i,j,k)/HH
                    k=k-1
                    if (k.ge.1) goto 111
                 end if
               end do
            end do
            do i=iimin-1,iimax
               do j=jjmin,jjmax
                  do k=1,kmax
                     huo(i,j,k)=min(hn(i,j,k),hn(i+1,j,k))
                     hun(i,j,k)=huo(i,j,k)
                  end do
                 kumin_pmz(i,j)=max(kmin_pmz(i,j),kmin_pmz(i+1,j))
               end do
            end do
            do i=iimin,iimax
               do j=jjmin-1,jjmax
                  do k=1,kmax
                     hvo(i,j,k)=min(hn(i,j,k),hn(i+1,j,k))
                     hvn(i,j,k)=hvo(i,j,k)
                  end do
                 kvmin_pmz(i,j)=max(kmin_pmz(i,j),kmin_pmz(i,j+1))
               end do
            end do
            do k=1,kmax
               do j=jjmin-1,jjmax
                  do i=iimin,iimax
                     hvo(i,j,k)=min(hn(i,j,k),hn(i+1,j,k))
                     hvn(i,j,k)=hvo(i,j,k)
                  end do
               end do
            end do
            kmin=1
            kumin=1
            kvmin=1
         case (3) ! general vertical coordinates
            do k=0,kmax
               ga(k) = k
            end do
            allocate(sig(0:kmax),stat=rc)    ! dimensionless sigma-coordinate
            if (rc /= 0) STOP 'coordinates: Error allocating (sig)'
            allocate(be(0:kmax),stat=rc)     ! dimensionless beta-coordinate
            if (rc /= 0) STOP 'coordinates: Error allocating (be)'
            allocate(gga(I3DFIELD),stat=rc)  ! dimensionless gamma-coordinate
            if (rc /= 0) stop 'coordinates: Error allocating memory (gga)'
            be(0)=  -_ONE_
            sig(0)= -_ONE_
            if (ddu .lt. _ZERO_) ddu=_ZERO_
            if (ddl .lt. _ZERO_) ddl=_ZERO_
            do k=1,kmax
               be(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
               be(k)=be(k)/(tanh(ddl)+tanh(ddu))-_ONE_
               sig(k)=k/float(kmax)-_ONE_
            end do
            if (gamma_surf) then
               kk=kmax
            else
               kk=1
            end if
            do j=jjmin-HALO,jjmax+HALO
               do i=iimin-HALO,iimax+HALO
                  HH=max(sseo(i,j)+H(i,j),min_depth)
                  alpha=min(&
                           ((be(kk)-be(kk-1))-D_gamma/HH&
                            *(sig(kk)-sig(kk-1)))&
                            /((be(kk)-be(kk-1))-(sig(kk)-sig(kk-1))),_ONE_)
                  gga(i,j,0)=-_ONE_
                  do k=1,kmax
                     gga(i,j,k)=alpha*sig(k)+(1.-alpha)*be(k)
                     if (gga(i,j,k) .lt. gga(i,j,k-1)) then
                        STDERR kk,(be(kk)-be(kk-1)),(sig(kk)-sig(kk-1))
                        STDERR D_gamma,HH
                        STDERR alpha
                        STDERR k-1,gga(i,j,k-1),be(k-1),sig(k-1)
                        STDERR k,gga(i,j,k),be(k),sig(k)
                        stop 'coordinates'
                     end if
                  end do
               end do
            end do

            kmin=1
            kumin=1
            kvmin=1
            kmin_pmz=1
            kumin_pmz=1
            kvmin_pmz=1

! Here, the initial layer distribution is calculated.
            do k=1,kmax
               do j=jjmin-HALO,jjmax+HALO
                  do i=iimin-HALO,iimax+HALO
                     HH=max(sseo(i,j)+H(i,j),min_depth)
                     hn(i,j,k)=HH*(gga(i,j,k)-gga(i,j,k-1))
                  end do
               end do
            end do

            do k=1,kmax
               do j=jjmin-HALO,jjmax+HALO
                  do i=iimin-HALO,iimax+HALO-1
                     HH=max(ssuo(i,j)+HU(i,j),min_depth)
                     huo(i,j,k)=HH*0.5*            &
                      (gga(i,j,k)-gga(i,j,k-1)+gga(i+1,j,k)-gga(i+1,j,k-1))
                     hun(i,j,k)=huo(i,j,k)
                  end do
               end do
            end do

!            allocate(zzz(I3DFIELD),stat=rc)  
!            do j=jjmin-HALO,jjmax+HALO
!               do i=iimin-HALO,iimax+HALO-1
!                  zzz(i,j,0)=-HU(i,j) 
!                  do k=1,kmax
!                     zzz(i,j,k)=zzz(i,j,k-1)+hun(i,j,k)
!                  end do
!               end do
!            end do
!            do j=jjmin-HALO,jjmax+HALO
!               do i=iimin-HALO,iimax+HALO
!                  do k=1,kmax
!                      write(90,*) (i-1)*0.7,zzz(i-1,j,k  )    
!                      write(90,*) (i-1)*0.7,zzz(i-1,j,k-1)    
!                      write(90,*) (i  )*0.7,zzz(i  ,j,k-1)    
!                      write(90,*) (i  )*0.7,zzz(i  ,j,k  )    
!                      write(90,*) (i-1)*0.7,zzz(i-1,j,k  )    
!                      write(90,*)     
!                      if ((zzz(i-1,j,k  ).lt.zzz(i  ,j,k-1)).or.(zzz(i-1,j,k-1).gt.zzz(i  ,j,k  )))                     &
!                 write(91,*) (i-0.5)*0.7,0.25*(zzz(i-1,j,k  )+ &
!                      zzz(i-1,j,k-1)+zzz(i  ,j,k-1)+zzz(i  ,j,k  )) 
!                  end do
!               end do
!            end do
!            stop

            do k=1,kmax
               do j=jjmin-HALO,jjmax+HALO-1
                  do i=iimin-HALO,iimax+HALO
                     HH=max(ssvo(i,j)+HV(i,j),min_depth)
                     hvo(i,j,k)=HH*0.5*            &
                      (gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
                     hvn(i,j,k)=hvo(i,j,k)
                  end do
               end do
            end do
!JMB
         case (4) ! sigma coordinates
         write(*,*) 'begin first '
            allocate(zpos(I3DFIELD),stat=rc)  ! z-coord. of interface
            if (rc /= 0) stop 'coordinates: Error allocating memory (zpos)'
            allocate(zposo(I3DFIELD),stat=rc)  ! z-coord. of interface
            if (rc /= 0) stop 'coordinates: Error allocating memory (zposo)'
            allocate(FF(I3DFIELD),stat=rc)  !
            if (rc /= 0) stop 'coordinates: Error allocating memory (FF)'
            allocate(work(I3DFIELD),stat=rc)  !
            if (rc /= 0) stop 'coordinates: Error allocating memory (work)'
            allocate(work2(I3DFIELD),stat=rc)  !
            if (rc /= 0) stop 'coordinates: Error allocating memory (work2)'
            allocate(work3(I3DFIELD),stat=rc)  !
            if (rc /= 0) stop 'coordinates: Error allocating memory (work3)'
            allocate(be(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (be)'
            allocate(jmb1(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (jmb1)'
            allocate(jmb2(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (jmb2)'
            allocate(jmb3(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (jmb3)'
            allocate(jmb4(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (jmb4)'
            allocate(NNloc(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (NNloc)'
            allocate(SSloc(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (SSloc)'
            allocate(av1(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (av1)'
            allocate(av2(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (av2)'
            allocate(av3(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (av3)'
            allocate(aav(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (aav)'
            allocate(gaa(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (gaa)'
            allocate(gaaold(0:kmax),stat=rc)     ! working space
            if (rc /= 0) STOP 'coordinates: Error allocating (gaaold)'
            kmaxm1= _ONE_/float(kmax)
! Dirty way to read initial distribution:

            do j=jjmin,jjmax
               do i=iimin,iimax
                  ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
                  hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
               end do
            end do
! Internal wave, only for testing for seiche test case
#ifdef 0           GOTO 10
            do j=jjmin,jjmax
                do i=iimin,iimax
                HH=H(i,j)*(0.5-0.25*(  &
                2*float(i)/float(iimax)-1 &
                ))
                DO K=1,kmax/2
                ii=kmax/2
                ho(i,j,k)=HH/float(ii)
                enddo
                do k=kmax/2,kmax
                ii=(kmax-kmax/2)
                ho(i,j,k)=(h(i,j)-HH)/float(ii)
                enddo
                enddo
            enddo
            hn=ho
 10       CONTINUE
#endif

            do j=jjmin,jjmax
               do i=iimin-1,iimax
                  huo(i,j,:)=(ssuo(i,j)+HU(i,j))*kmaxm1
                  hun(i,j,:)=(ssun(i,j)+HU(i,j))*kmaxm1
               end do
            end do
            do j=jjmin-1,jjmax
               do i=iimin,iimax
                  hvo(i,j,:)=(ssvo(i,j)+HV(i,j))*kmaxm1
                  hvn(i,j,:)=(ssvn(i,j)+HV(i,j))*kmaxm1
               end do
            end do
            kmin=1
            kumin=1
            kvmin=1

!END JMB


         case default

      end select

   end if ! first

   select case (cord_type)
      case (1) ! sigma coordinates
         if (equiv_sigma) then
            kmaxm1= _ONE_/float(kmax)
            do j=jjmin-HALO,jjmax+HALO
               do i=iimin-HALO,iimax+HALO
                  ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
                  hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
               end do
            end do

            do j=jjmin-HALO,jjmax+HALO
               do i=iimin-HALO,iimax+HALO-1
                  huo(i,j,:)=(ssuo(i,j)+HU(i,j))*kmaxm1
                  hun(i,j,:)=(ssun(i,j)+HU(i,j))*kmaxm1
               end do
            end do

            do j=jjmin-HALO,jjmax+HALO-1
               do i=iimin-HALO,iimax+HALO
                  hvo(i,j,:)=(ssvo(i,j)+HV(i,j))*kmaxm1
                  hvn(i,j,:)=(ssvn(i,j)+HV(i,j))*kmaxm1
               end do
            end do

         else ! non-equivdistant

            do j=jjmin-HALO,jjmax+HALO
               do i=iimin-HALO,iimax+HALO
                  ho(i,j,1:kmax)=(sseo(i,j)+H(i,j))*dga(1:kmax)
                  hn(i,j,1:kmax)=(ssen(i,j)+H(i,j))*dga(1:kmax)
               end do
            end do

            do j=jjmin-HALO,jjmax+HALO
               do i=iimin-HALO,iimax+HALO-1
                  huo(i,j,1:kmax)=(ssuo(i,j)+HU(i,j))*dga(1:kmax)
                  hun(i,j,1:kmax)=(ssun(i,j)+HU(i,j))*dga(1:kmax)
               end do
            end do

            do j=jjmin-HALO,jjmax+HALO-1
               do i=iimin-HALO,iimax+HALO
                  hvo(i,j,1:kmax)=(ssvo(i,j)+HV(i,j))*dga(1:kmax)
                  hvn(i,j,1:kmax)=(ssvn(i,j)+HV(i,j))*dga(1:kmax)
               end do
            end do
         end if

      case (2) ! z-level
!        Maybe use mask and not loop over k, e.g. ho(i,j,:)=hn(i,j,:)
!kbk#define USE_ARRAY_SYNTAX
         do j=jjmin-HALO,jjmax+HALO
            do i=iimin-HALO,iimax+HALO
               HH=ssen(i,j)+H(i,j)
#ifdef USE_ARRAY_SYNTAX
               ho(i,j,:)=hn(i,j,:)
               hn(i,j,:)=gga(i,j,:)*HH
#else
               do k=1,kmax
                  ho(i,j,k)=hn(i,j,k)
                  hn(i,j,k)=gga(i,j,k)*HH
               end do
#endif
            end do
         end do
         do j=jjmin-HALO,jjmax+HALO
            do i=iimin-HALO,iimax+HALO-1
               HH=ssun(i,j)+HU(i,j)
#ifdef USE_ARRAY_SYNTAX
               huo(i,j,:)=hun(i,j,:)
               hun(i,j,:)=0.5*(gga(i,j,:)+gga(i+1,j,:))*HH
#else
               do k=1,kmax
                  huo(i,j,k)=hun(i,j,k)
                  hun(i,j,k)=0.5*(gga(i,j,k)+gga(i+1,j,k))*HH
               end do
#endif
            end do
         end do
         do j=jjmin-HALO,jjmax+HALO-1
            do i=iimin-HALO,iimax+HALO
               HH=ssvn(i,j)+HV(i,j)
#ifdef USE_ARRAY_SYNTAX
               hvo(i,j,:)=hvn(i,j,:)
               hvn(i,j,:)=0.5*(gga(i,j,:)+gga(i,j+1,:))*HH
#else
               do k=1,kmax
                  hvo(i,j,k)=hvn(i,j,k)
                  hvn(i,j,k)=0.5*(gga(i,j,k)+gga(i,j+1,k))*HH
               end do
#endif

#undef USE_ARRAY_SYNTAX
            end do
         end do

      case (3) ! general vertical coordinates

! The general vertical coordinates can be relaxed towards the new layer
! thicknesses by the following relaxation time scale r. This should
! later be generalised also for sigma coordinates.

         do j=jjmin-HALO,jjmax+HALO
            do i=iimin-HALO,iimax+HALO
               r=cord_relax/dt*H(i,j)/maxdepth
               HH=ssen(i,j)+H(i,j)
               if (HH .lt. D_gamma) then
                  do k=1,kmax
                     ho(i,j,k)=hn(i,j,k)
                     hn(i,j,k)=HH/float(kmax)
                  end do
               else
                  zz=-H(i,j)
                  do k=1,kmax-1
                     ho(i,j,k)=hn(i,j,k)
                     hn(i,j,k)=(ho(i,j,k)*r+HH*(gga(i,j,k)-gga(i,j,k-1)))/(r+1.)
                     zz=zz+hn(i,j,k)
                  end do
                  ho(i,j,kmax)=hn(i,j,kmax)
                  hn(i,j,kmax)=ssen(i,j)-zz
               end if
            end do
         end do

         do j=jjmin-HALO,jjmax+HALO
            do i=iimin-HALO,iimax+HALO-1
!KBK               if (au(i,j) .gt. 0) then
                  r=cord_relax/dt*HU(i,j)/maxdepth
                  zz=-HU(i,j)
                  HH=ssun(i,j)+HU(i,j)
                  do k=1,kmax-1
                     huo(i,j,k)=hun(i,j,k)
                     hun(i,j,k)=(huo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                               +gga(i+1,j,k)-gga(i+1,j,k-1)))/(r+1.)
                     zz=zz+hun(i,j,k)
                  end do
                  huo(i,j,kmax)=hun(i,j,kmax)
                  hun(i,j,kmax)=ssun(i,j)-zz
!KBK               end if
            end do
         end do

         do j=jjmin-HALO,jjmax+HALO-1
            do i=iimin-HALO,iimax+HALO
!KBK               if (av(i,j).gt.0) then
                  r=cord_relax/dt*HV(i,j)/maxdepth
                  zz=-HV(i,j)
                  HH=ssvn(i,j)+HV(i,j)
                  do k=1,kmax-1
                     hvo(i,j,k)=hvn(i,j,k)
                     hvn(i,j,k)=(hvo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                               +gga(i,j+1,k)-gga(i,j+1,k-1)))/(r+1.)
                     zz=zz+hvn(i,j,k)
                  end do
                  hvo(i,j,kmax)=hvn(i,j,kmax)
                  hvn(i,j,kmax)=ssvn(i,j)-zz
!KBK               end if
            end do
         end do
!JMB
! DX and DY scaling should be included....
      case (4) ! adaptive grids
         write(*,*) 'begin adaptive grids '
         ho=hn


! Lagrangien step 
           do k=1,kmax
              do j=jjmin+1,jjmax-1
                 do i=iimin+1,iimax-1
                    hn(i,j,k)=ho(i,j,k)                      &
                        -((uu(i,j,k)*DYU-uu(i-1,j  ,k)*DYUIM1) &
                        +(vv(i,j,k)*DXV-vv(i  ,j-1,k)*DXVJM1))*ARCD1*dt*FACLAG &
                        +( &
                           (ho(i+1,j,k)-ho(i,j,k))*au(i,j) &
                          +(ho(i-1,j,k)-ho(i,j,k))*au(i-1,j) &
                          +(ho(i,j+1,k)-ho(i,j,k))*av(i,j) &
                          +(ho(i,j-1,k)-ho(i,j,k))*av(i,j-1) &
                         )*FACDIF*0.25
! Ensure smooth transition to cut-off around depth
                    hn(i,j,k)=max(hn(i,j,k),depthmin/5.)
                    IF(hn(i,j,k).LT.depthmin) hn(i,j,k)=hn(i,j,k)+depthmin/5.
                 end do
              end do
           end do

           call hcheck(hn,ssen,h)
! First step provided Lagrangian modified hn with additional horizontal
! diffusion of hn. Consistency with total depth is ensured




! Horizontal diffusion of zpos repeated MHOR times
         do ii=1,MHOR
! Prepare zpos field
         call htoz(hn,zpos)

! Next line probably not needed
!        call htoz(ho,zposo)


! For problem zones, put zero 'horizontal diffusion'
            do k=1,kmax
               do j=jjmin,jjmax
                  do i=iimin,iimax
                     work2(i,j,k)=FACHOR*au(i,j)
                  enddo
                enddo
             enddo
             do k=2,kmax
                do j=jjmin,jjmax
                   do i=iimin+1,iimax
                      if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                         work2(i,j,k)=0
                         work2(i,j,k-1)=0
                         work2(i-1,j,k)=0
                         work2(i-1,j,k-1)=0
                      endif
                    enddo
                 enddo
              enddo
              do   k=1,kmax
                 do j=jjmin,jjmax
                    do i=iimin,iimax
                       work3(i,j,k)=FACHOR*av(i,j)
                     enddo
                 enddo
              enddo
              do k=2,kmax
                 do j=jjmin+1,jjmax
                    do i=iimin,iimax
                       if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                          work3(i,j,k)=0
                          work3(i,j,k-1)=0
                          work3(i,j-1,k)=0
                          work3(i,j-1,k-1)=0
                       endif
                    enddo
                 enddo
              enddo
! Check for  bounds on i j etc and possible dirichlet conditions on z

! Dirty BC
           work=zpos
             do k=1,kmax-1
                do j=jjmin+1,jjmax-1
                   do i=iimin+1,iimax-1
                   RM=0
                   IM=0

                      DO iii=max(iimin,i-iw),min(iimax,i+iw)
                        DO jjj=max(jjmin,j-iw),min(jjmax,j+iw)
                           RM=RM+az(iii,jjj)*(rho(III,jjj,k+1)+rho(III,jjj,k))
                           IM=IM+az(iii,jjj)
                        ENDDO
                      ENDDO

                      work(i,j,k)=zpos(i,j,k)+ &
                              (   &
                           (zpos(i+1,j,k)-zpos(i,j,k))*work2(i,j,k) - &
                           (zpos(i,j,k)-zpos(i-1,j,k))*work2(i-1,j,k)  &
                         + (zpos(i,j+1,k)-zpos(i,j,k))*work3(i,j,k) - &
                           (zpos(i,j,k)-zpos(i,j-1,k))*work3(i,j-1,k)  &
                              )*0.25    &
                             +( &
                             .5* (RM/IM-rho(i,j,k+1)-rho(i,j,k))/NN(i,j,k)&
                              )*FACISO



                   end do
                end do
             end do
! Local consistency
            do k=1,kmax
               do j=jjmin,jjmax
                  do i=iimin,iimax
                     zpos(i,j,k)=work(i,j,k)
                     if((zpos(i,j,k)-zpos(i,j,k-1)).lt.depthmin) then
                        zpos(i,j,k)=zpos(i,j,k-1)+depthmin
                     endif
                  enddo
               enddo
            enddo
 
            call ztoh(zpos,hn,depthmin)
            call hcheck(hn,ssen,h)


         end do ! End of Horizontal diffusion of zpos repeated MHOR times
! After horizontal diffusion updated new and depth consistent hn is available


#if 0
     do ii=1,MVER
! Interpolate to new positions
         call htoz(hn,zpos)
         call htoz(ho,zposo)
          do j=jjmin,jjmax
             do i=iimin,iimax
               do k=1,kmax
                  jmb1(k)=(zpos(i,j,k)+zpos(i,j,k-1))*0.5
                  jmb2(k)=(zposo(i,j,k)+zposo(i,j,k-1))*0.5
                  jmb3(k)=S(i,j,k)
                enddo
!               call VINTEP(jmb3,jmb2,0.,0.,1,KMAX,jmb1,1,KMAX,jmb4)
! col_interpol(N,cols,obs_z,obs_prof,nlev,model_z,model_prof)
                call col_interpol(kmax,1,jmb2,jmb3,kmax,jmb1,jmb4)

!                if ((i.eq.10).and.(j.eq.2))
!               do k=1,kmax
!                write(*,*) k,zposo(i,j,k),jmb   
!                enddo
!                end if
                do k=1,kmax
                   FF(i,j,k)=jmb4(k)
!                  FF(i,j,k)=S(i,j,k)
                enddo
             enddo
          enddo

          work=zpos
          jmb1=0
          do j=jjmin,jjmax
             do i=iimin,iimax
                do k=1,kmax-1
                   jmb1(k)=min(ABS(FF(i,j,k+1)-FF(i,j,k))/(zpos(i,j,k+1)-zpos(i,j,k)),30.)+0.00000

! TEST at the left boundary
               enddo
               jmb1(0)=jmb1(1)
               jmb1(kmax)=jmb1(kmax-1)
! One iteration does not need reinterpolation...
               do k=1,kmax-1
                  zpos(i,j,k)=work(i,j,k)+ FACVER*(  &
                     (work(i,j,k+1)-work(i,j,k))*(jmb1(k+1)+jmb1(k)) &
                     -(work(i,j,k)-work(i,j,k-1))*(jmb1(k-1)+jmb1(k)) &
                                    )
               enddo
             enddo
           enddo

! End vertical loop
      enddo
! End Test of vertical diffusion of zpos

#else
         call htoz(hn,zpos)
         dtgrid=dt/float(split)
          do j=jjmin,jjmax
             do i=iimin,iimax
                NNloc=NN(i,j,:) 
                SSloc=SS(i,j,:) 
                do k=0,kmax
                   gaa(k)=(zpos(i,j,k)-ssen(i,j))/(ssen(i,j)+H(i,j))
                   gaaold(k)=gaa(k)
                end do 
                do ii=1,split
!     Stratification
! Dirty jm
      NNloc(kmax)=NNloc(kmax-1)       
      NNloc(0)=NNloc(1)       
      SSloc(kmax)=SSloc(kmax-1)       
      SSloc(0)=SSloc(1)       
      do k=1,kmax
         av1(k)=min(1.,max(0.,0.5*(NNloc(k)+NNloc(k-1)))/g*1024./NNnorm)
      end do

!     Shear
      do k=1,kmax
         av2(k)=min(1.,sqrt(max(0.,0.5*(SSloc(k)+SSloc(k-1))))/SSnorm)
      end do

!     Distance from surface
      do k=1,kmax
         av3(k)=1./(dsurf-0.5*(gaa(k-1)+gaa(k))*H(i,j))+            &
                1./(dsurf+(1.+0.5*(gaa(k-1)+gaa(k)))*H(i,j))
      end do

!     Calculation of grid diffusivity
      do k=1,kmax
         aav(k)=H(i,j)/Tgrid*(c1ad*av1(k)+c2ad*av2(k)+c3ad*av3(k)+c4ad/H(i,j))
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
             enddo
          enddo
#endif


! Back to hn
      call ztoh(zpos,hn,depthmin)
! Normally if positive defined vertical diffusion no check required


! Finally derive interface grid sizes for uu and vv
! Interface treatment and check
! uu
            huo=hun
            do k=1,kmax
               do j=jjmin,jjmax
                  do i=iimin,iimax-1
                     hun(i,j,k)=0.5*(hn(i,j,k)+hn(i+1,j,k))
                  end do
               end do
            end do
! maybe not allowed in iimax....
            call hcheck(hun,ssun,hu)
! vv
            hvo=hvn
            do k=1,kmax
               do j=jjmin,jjmax-1
                  do i=iimin,iimax
                     hvn(i,j,k)=0.5*(hn(i,j,k)+hn(i,j+1,k))
                  end do
               end do
            end do
! maybe not allowed in jjmax....
            call hcheck(hvn,ssvn,hv)



      case default

   end select

#ifdef DEBUG
   write(debug,*) 'Leaving coordinates()'
   write(debug,*)
#endif
   return
   end subroutine coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
