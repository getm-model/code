!$Id: slow_diffusion.F90,v 1.4 2003-08-28 15:19:03 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_diffusion() - diffusion of momentum.
!
! !INTERFACE:
   subroutine slow_diffusion(AM)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,az,au,av,ax,H,HU,HV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: D,U,V,UEx,VEx,Uint,Vint,PP
   use variables_3d, only: ssen,ssun,ssvn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  REALTYPE, intent(in)                 :: AM
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii,jj
   REALTYPE                  :: Di(iimin-1:iimax+1,jjmin-1:jjmax+1)
   REALTYPE                  :: DUi(iimin-1:iimax+1,jjmin-1:jjmax+1)
   REALTYPE                  :: DVi(iimin-1:iimax+1,jjmin-1:jjmax+1)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_diffusion() # ',Ncall
#endif

   do j=jjmin-1,jjmax+1
      do i=iimin-1,iimax+1
         Di(i,j)=ssen(i,j)+H(i,j)
      end do
   end do

   do j=jjmin-1,jjmax+1
      do i=iimin-1,iimax+1
         DUi(i,j)=ssun(i,j)+HU(i,j)
      end do
   end do

   do j=jjmin-1,jjmax+1
      do i=iimin-1,iimax+1
         DVi(i,j)=ssvn(i,j)+HV(i,j)
      end do
   end do

! Central for dx(2*AM*dx(U^2/HU))
   do j=jjmin,jjmax
      do i=iimin,iimax+1          ! PP defined on T-points
         if (az(i,j) .ge. 1) then
            PP(i,j)=2.*AM*DYC*Di(i,j)               &
               *(Uint(i,j)/DUi(i,j)-Uint(i-1,j)/DUi(i-1,j))/DXC
         else
            PP(i,j)=_ZERO_
         end if
      end do
   end do
   do j=jjmin,jjmax      ! UEx defined on U-points
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=UEx(i,j)-(PP(i+1,j)-PP(i  ,j))*ARUD1
         end if
      end do
   end do

! Central for dy(AM*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jjmin-1,jjmax        ! PP defined on X-points
      do i=iimin,iimax
         if (ax(i,j) .ge. 1) then
            PP(i,j)=AM*0.5*(DUi(i,j)+DUi(i,j+1))*DXX  &
                   *((Uint(i,j+1)/DUi(i,j+1)-Uint(i,j)/DUi(i,j))/DYX &
                    +(Vint(i+1,j)/DVi(i+1,j)-Vint(i,j)/DVi(i,j))/DXX )
         else
            PP(i,j)=_ZERO_
         end if
      end do
   end do
   do j=jjmin,jjmax        !UEx defined on U-points
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=UEx(i,j)-(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do

! Central for dx(AM*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jjmin,jjmax      ! PP defined on X-points
      do i=iimin-1,iimax
         if (ax(i,j) .ge. 1) then
            PP(i,j)=AM*0.5*(DVi(i,j)+DVi(i+1,j))*DXX  &
                   *((Uint(i,j+1)/DUi(i,j+1)-Uint(i,j)/DUi(i,j))/DYX &
                    +(Vint(i+1,j)/DVi(i+1,j)-Vint(i,j)/DVi(i,j))/DXX )
         else
            PP(i,j)=_ZERO_
         end if
      end do
   end do
   do j=jjmin,jjmax          ! VEx defined on V-points
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=VEx(i,j)-(PP(i  ,j)-PP(i-1,j))*ARVD1
         end if
      end do
   end do

! Central for dy(2*AM*dy(V^2/DV))
   do j=jjmin,jjmax+1     ! PP defined on T-points
      do i=iimin,iimax
         if (az(i,j) .ge. 1) then
            PP(i,j)=2.*AM*DXC*Di(i,j)               &
                   *(Vint(i,j)/DVi(i,j)-Vint(i,j-1)/DVi(i,j-1))/DYC
         else
            PP(i,j)=_ZERO_
         end if
      end do
   end do
   do j=jjmin,jjmax             ! VEx defined on V-points
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=VEx(i,j)-(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do

#ifdef DEBUG
     write(debug,*) 'Leaving slow_diffusion()'
     write(debug,*)
#endif
   return
   end subroutine slow_diffusion
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
