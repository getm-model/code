!$Id: ip_blumberg_mellor.F90,v 1.2 2004-04-20 15:52:58 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_blumberg_mellor()
!
! !INTERFACE:
   subroutine ip_blumberg_mellor()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,au,av,H,HU,HV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyv
#else
   use domain, only: dx,dy
#endif
   use variables_3d, only: kmin,hn,hun,hvn,idpdx,idpdy,rho
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: ip_blumberg_mellor.F90,v $
!  Revision 1.2  2004-04-20 15:52:58  hb
!  missing factor 0.5 added
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/01 15:50:13  gotm
!  removed dead print statement
!
!  Revision 1.1.1.1  2002/05/02 14:00:59  gotm
!  recovering after CVS crash
!
!  Revision 1.15  2001/10/22 08:23:53  bbh
!  Removed reference to kplus and kminus when not -DPRESS_GRAD_Z
!
!  Revision 1.14  2001/10/22 07:47:59  bbh
!  Needed to run dos2unix
!
!  Revision 1.13  2001/10/22 07:45:26  bbh
!  Fixed a serious bug - now gives the same as old version
!
!  Revision 1.10  2001/09/04 07:29:18  bbh
!  Internal pressure based on interpolation to z-coordinates
!
!  Revision 1.9  2001/09/03 13:03:37  bbh
!  Initial pressure gradient can now be subtracted
!
!  Revision 1.8  2001/08/31 15:40:37  bbh
!  initial pressure can be subtracted now
!
!  Revision 1.7  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/25 19:09:48  bbh
!  Removed dead code
!
!  Revision 1.4  2001/05/20 09:19:09  bbh
!  Use who specified twice
!
!  Revision 1.3  2001/05/20 07:51:40  bbh
!  Internal pressure included
!
!  Revision 1.2  2001/05/11 13:47:00  bbh
!  Added actual code
!

!  Added further support for baroclinicity
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,rc
   REALTYPE                  :: dxm1,dym1,x,y,x1,y1,hc
   REALTYPE                  :: grdl,grdu,rhol,rhou,prgr,dxz,dyz
   LOGICAL,save              :: first=.true.
   REALTYPE,dimension(:,:,:), allocatable        :: zz
   REALTYPE,save,dimension(:,:,:), allocatable   :: idpdx0
   REALTYPE,save,dimension(:,:,:), allocatable   :: idpdy0
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ip_blumberg_mellor() # ',Ncall
#endif

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

#ifdef SUBSTR_INI_PRESS
   if (first) then
      allocate(idpdx0(I3DFIELD),stat=rc)    ! Initial x - pressure gradient.
      if (rc /= 0) stop 'ip_blumberg_mellor.F90: Error allocating memory (idpdx0)'
      allocate(idpdy0(I3DFIELD),stat=rc)    ! Initial y - pressure gradient.
      if (rc /= 0) stop 'ip_blumberg_mellor.F90: Error allocating memory (idpdy0)'
   end if
#endif

   allocate(zz(I3DFIELD),stat=rc)    ! Interface heights.
   if (rc /= 0) stop 'ip_blumberg_mellor.F90: Error allocating memory (zz)'

   zz = _ZERO_
   idpdx = _ZERO_
   idpdy = _ZERO_

! First, the interface heights are calculated in order to get the
! interface slopes further down.
   do j=jjmin,jjmax+1
      do i=iimin,iimax+1
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+hn(i,j,1)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+hn(i,j,k)
            end do
         end if
      end do
   end do

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the u-velocity equation.
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            grdl=(rho(i+1,j,kmax)-rho(i,j,kmax))*dxm1
            rhol=0.5*(rho(i+1,j,kmax)+rho(i,j,kmax))
            prgr=grdl
            idpdx(i,j,kmax)=hun(i,j,kmax)*prgr*0.5*hun(i,j,kmax)
            do k=kmax-1,1,-1
               grdu=grdl
               grdl=(rho(i+1,j,k)-rho(i,j,k))*dxm1
               rhou=rhol
               rhol=0.5*(rho(i+1,j,k)+rho(i,j,k))
               dxz=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               prgr=prgr+0.5*(hun(i,j,k)+hun(i,j,k+1))*0.5*(grdu+grdl)-dxz*(rhou-rhol)
!               hc=(zz(i+1,j,k)-zz(i+1,j,k+1))**2*(rho(i+1,j,k)-rho(i+1,j,k+1))
!               hc=hc-(zz(i,j,k)-zz(i,j,k+1))**2*(rho(i,j,k)-rho(i,j,k+1)) 
!               hc=hc+(zz(i+1,j,k+1)-zz(i,j,k+1))**2*(rho(i+1,j,k+1)-rho(i,j,k+1)) 
!               hc=hc-(zz(i+1,j,k)-zz(i,j,k))**2*(rho(i+1,j,k)-rho(i,j,k)) 
!               hc=hc*dxm1*0.1666666667/(zz(i,j,k)+zz(i+1,j,k)-zz(i,j,k+1)-zz(i+1,j,k+1))
!               prgr=prgr+0.5*(hun(i,j,k)+hun(i,j,k+1))*hc
               idpdx(i,j,k)=hun(i,j,k)*prgr
            end do
         end if
      end do
   end do

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dym1 = _ONE_/DYV
#endif
            grdl=(rho(i,j+1,kmax)-rho(i,j,kmax))*dym1
            rhol=0.5*(rho(i,j+1,kmax)+rho(i,j,kmax))
            prgr=grdl
            idpdy(i,j,kmax)=hvn(i,j,kmax)*prgr*0.5*hvn(i,j,kmax)
            do k=kmax-1,1,-1
               grdu=grdl
               grdl=(rho(i,j+1,k)-rho(i,j,k))*dym1
               rhou=rhol
               rhol=0.5*(rho(i,j+1,k)+rho(i,j,k))
               dyz=(zz(i,j+1,k)-zz(i,j,k))*dym1
               prgr=prgr+0.5*(hvn(i,j,k)+hvn(i,j,k+1))*0.5*(grdu+grdl)-dyz*(rhou-rhol)
!               hc=(zz(i,j+1,k)-zz(i,j+1,k+1))**2*(rho(i,j+1,k)-rho(i,j+1,k+1))
!               hc=hc-(zz(i,j,k)-zz(i,j,k+1))**2*(rho(i,j,k)-rho(i,j,k+1)) 
!               hc=hc+(zz(i,j+1,k+1)-zz(i,j,k+1))**2*(rho(i,j+1,k+1)-rho(i,j,k+1)) 
!               hc=hc-(zz(i,j+1,k)-zz(i,j,k))**2*(rho(i,j+1,k)-rho(i,j,k)) 
!               hc=hc*dym1*0.1666666667/(zz(i,j,k)+zz(i,j+1,k)-zz(i,j,k+1)-zz(i,j+1,k+1))
!               prgr=prgr+0.5*(hvn(i,j,k)+hvn(i,j,k+1))*hc
               idpdy(i,j,k)=hvn(i,j,k)*prgr
            end do
         end if
      end do
   end do

#ifdef SUBSTR_INI_PRESS
   if (first) then
      first = .false.
      do k=0,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
               idpdx0(i,j,k) = idpdx(i,j,k)
               idpdx(i,j,k) = _ZERO_
               idpdy0(i,j,k) = idpdy(i,j,k)
               idpdy(i,j,k) = _ZERO_
            end do
         end do
      end do
   else
      do k=0,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
               idpdx(i,j,k) = idpdx(i,j,k) - idpdx0(i,j,k)
               idpdy(i,j,k) = idpdy(i,j,k) - idpdy0(i,j,k)
            end do
         end do
      end do
   end if
#endif

#ifdef FORTRAN90
   deallocate(zz,stat=rc)    ! work array
   if (rc /= 0) stop 'ip_blumberg_mellor: Error de-allocating memory (zz)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving ip_blumberg_mellor()'
   write(debug,*)
#endif
   return
   end subroutine ip_blumberg_mellor
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
