!$Id: coordinates.F90,v 1.3 2003-04-23 12:16:27 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  coordinates() - defines the vertical coordinate system.
!
! !INTERFACE:
   subroutine coordinates(cord_type,cord_relax)
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG,HU_TAG,HV_TAG
   use domain,   only: iimin,iimax,jjmin,jjmax,kmax,H,HU,HV,az,au,av,min_depth
   use domain,   only: ga,ddu,ddl,d_gamma
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hun,hvo,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: cord_type
   REALTYPE, intent(in)                :: cord_relax
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: coordinates.F90,v $
!  Revision 1.3  2003-04-23 12:16:27  kbk
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
   integer :: i,j,k,rc,kk
   REALTYPE:: tmp,kmaxm1,alpha,maxdepth
   LOGICAL   :: gamma_surf=.true.
   REALTYPE:: HH,zz,r
   logical, save :: first=.true.,equiv_sigma=.false.
   REALTYPE, save, dimension(:),     allocatable  :: dga,be,sig
   REALTYPE, save, dimension(:,:,:), allocatable  :: gga
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
      allocate(ga(0:kmax),stat=rc)
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
         case (2) ! z-level
            STDERR 'z-level not coded yet'
            stop 'coordinates'
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
            do j=jjmin,jjmax
               do i=iimin,iimax
                  HH=max(H(i,j),min_depth)
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

            gga(iimin:iimax,jjmin-1,:) = gga(iimin:iimax,jjmin,:)
            gga(iimin:iimax,jjmax+1,:) = gga(iimin:iimax,jjmax,:)
            gga(iimin-1,jjmin:jjmax,:)=gga(iimin,jjmin:jjmax,:)
            gga(iimax+1,jjmin:jjmax,:)=gga(iimax,jjmin:jjmax,:)

            call update_3d_halo(gga,gga,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
            call wait_halo(H_TAG)
            kmin=1
            kumin=1
            kvmin=1

! Here, the initial layer distribution is calculated.

         do k=1,kmax
            do j=jjmin,jjmax
               do i=iimin,iimax
                  HH=max(sseo(i,j)+H(i,j),min_depth)
                  hn(i,j,k)=HH*(gga(i,j,k)-gga(i,j,k-1))
               end do
            end do
         end do
         do k=1,kmax
            do j=jjmin,jjmax
               do i=iimin-1,iimax
                  HH=max(ssuo(i,j)+HU(i,j),min_depth)
                  huo(i,j,k)=HH*0.5*            &
                   (gga(i,j,k)-gga(i,j,k-1)+gga(i+1,j,k)-gga(i+1,j,k-1))
                  hun(i,j,k)=huo(i,j,k)
               end do
            end do
         end do

         do k=1,kmax
            do j=jjmin-1,jjmax
               do i=iimin,iimax
                  HH=max(ssvo(i,j)+HV(i,j),min_depth)
                  hvo(i,j,k)=HH*0.5*            &
                   (gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
                  hvn(i,j,k)=hvo(i,j,k)
               end do
            end do
         end do


         case default
      end select

   end if ! first

   select case (cord_type)
      case (1) ! sigma coordinates
         if (equiv_sigma) then
            kmaxm1= _ONE_/float(kmax)
            do j=jjmin,jjmax
               do i=iimin,iimax
                  ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
                  hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
               end do
            end do

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

         else ! non-equivdistant
            do j=jjmin,jjmax
               do i=iimin,iimax
                  ho(i,j,1:kmax)=(sseo(i,j)+H(i,j))*dga(1:kmax)
                  hn(i,j,1:kmax)=(ssen(i,j)+H(i,j))*dga(1:kmax)
               end do
            end do

            do j=jjmin,jjmax
               do i=iimin-1,iimax
                  huo(i,j,1:kmax)=(ssuo(i,j)+HU(i,j))*dga(1:kmax)
                  hun(i,j,1:kmax)=(ssun(i,j)+HU(i,j))*dga(1:kmax)
               end do
            end do

            do j=jjmin-1,jjmax
               do i=iimin,iimax
                  hvo(i,j,1:kmax)=(ssvo(i,j)+HV(i,j))*dga(1:kmax)
                  hvn(i,j,1:kmax)=(ssvn(i,j)+HV(i,j))*dga(1:kmax)
               end do
            end do
         end if
      case (2) ! z-level
      case (3) ! general vertical coordinates

! The general vertical coordinates can be relaxed towards the new layer
! thicknesses by the following relaxation time scale r. This should
! later be generalised also for sigma coordinates.

         maxdepth=600. ! needs to be calculated later ... HB 
         do j=jjmin,jjmax
            do i=iimin,iimax
               r = cord_relax/dt*H(i,j)/maxdepth
               zz = -H(i,j)
               do k=1,kmax-1
                  ho(i,j,k) = hn(i,j,k)
                  HH=max(ssen(i,j)+H(i,j),min_depth)
                  hn(i,j,k)=(ho(i,j,k)*r+HH*(gga(i,j,k)-gga(i,j,k-1)))/(r+1.)
                  zz = zz + hn(i,j,k)
               end do
               ho(i,j,kmax) = hn(i,j,kmax)
               hn(i,j,kmax)=ssen(i,j)-zz
            end do
         end do

         do j=jjmin,jjmax
            do i=iimin-1,iimax
               if (au(i,j) .gt. 0) then
                  r = cord_relax/dt*HU(i,j)/maxdepth
                  zz = -HU(i,j)
                  do k=1,kmax-1
                     huo(i,j,k) = hun(i,j,k)
                     HH=max(ssun(i,j)+HU(i,j),min_depth)
                     hun(i,j,k)=(huo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                               +gga(i+1,j,k)-gga(i+1,j,k-1)))/(r+1.)
                     zz = zz + hun(i,j,k)
                  end do
                  huo(i,j,kmax) = hun(i,j,kmax)
                  hun(i,j,kmax)=ssun(i,j)-zz
               end if
            end do
         end do

         do j=jjmin-1,jjmax
            do i=iimin,iimax
               if (av(i,j).gt.0) then
                  r = cord_relax/dt*HV(i,j)/maxdepth
                  zz = -HV(i,j)
                  do k=1,kmax-1
                     hvo(i,j,k) = hvn(i,j,k)
                     HH=max(ssvn(i,j)+HV(i,j),min_depth)
                     hvn(i,j,k)=(hvo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                               +gga(i,j+1,k)-gga(i,j+1,k-1)))/(r+1.)
                     zz=zz+hvn(i,j,k)
                  end do
                  hvo(i,j,kmax) = hvn(i,j,kmax)
                  hvn(i,j,kmax)=ssvn(i,j)-zz
               end if
            end do
         end do

      case default

   end select
   call update_3d_halo(ho,ho,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
   call update_3d_halo(hn,hn,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
   call update_3d_halo(huo,huo,au,iimin,jjmin,iimax,jjmax,kmax,HU_TAG)
   call wait_halo(HU_TAG)
   call update_3d_halo(hun,hun,au,iimin,jjmin,iimax,jjmax,kmax,HU_TAG)
   call wait_halo(HU_TAG)
   call update_3d_halo(hvo,hvo,av,iimin,jjmin,iimax,jjmax,kmax,HV_TAG)
   call wait_halo(HV_TAG)
   call update_3d_halo(hvn,hvn,av,iimin,jjmin,iimax,jjmax,kmax,HV_TAG)
   call wait_halo(HV_TAG)


#ifdef DEBUG
   write(debug,*) 'Leaving Coordinates()'
   write(debug,*)
#endif
   return
   end subroutine coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
