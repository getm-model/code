!$Id: ip_blumberg_mellor_lin.F90,v 1.1 2004-04-06 12:42:50 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_blumberg_mellor_lin()
!
! !INTERFACE:
   subroutine ip_blumberg_mellor_lin()
!
! !DESCRIPTION:
!
! !USES:
   use internal_pressure
   use variables_3d, only: kumin_pmz,kvmin_pmz
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: ip_blumberg_mellor_lin.F90,v $
!  Revision 1.1  2004-04-06 12:42:50  kbk
!  internal pressure calculations now uses wrapper
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1
   REALTYPE                  :: prgr,dxzu,dxzl,dyzu,dyzl
   REALTYPE                  :: dzr2,dzr1,dxru,dxrl,dyru,dyrl,aa,bb,cc
!EOP
!-----------------------------------------------------------------------
!BOC
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the pressure point heights are calculated in order to get the
!  interface slopes further down.
   do j=jjmin,jjmax+1
      do i=iimin,iimax+1
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+0.5*hn(i,j,1)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+0.5*(hn(i,j,k-1)+hn(i,j,k))
            end do
         end if
      end do
   end do

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the u-velocity equation.
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            dxzl=(zz(i+1,j,kmax)-zz(i,j,kmax))*dxm1
            dxrl=(rho(i+1,j,kmax)-rho(i,j,kmax))*dxm1
            prgr=0.
            do k=kmax-1,kumin_pmz(i,j),-1
               dxzu=dxzl
               dxzl=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               dxru=dxrl
               dxrl=(rho(i+1,j,k)-rho(i,j,k))*dxm1
               dzr2=(rho(i+1,j,k+1)-rho(i+1,j,k))/(zz(i+1,j,k+1)-zz(i+1,j,k))
               dzr1=(rho(i  ,j,k+1)-rho(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
               aa=0.5*(dxrl+dxru)
               bb=0.5*(dxzl+dxzu)
               cc=0.5*(dzr2+dzr1)
               prgr=prgr+(aa-bb*cc)*0.5*(hun(i,j,k+1)+hun(i,j,k))
               if (k.eq.kmax-1) idpdx(i,j,kmax)=prgr*0.5*hun(i,j,kmax)**2
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
           dyzl=(zz(i,j+1,kmax)-zz(i,j,kmax))*dym1
           dyrl=(rho(i,j+1,kmax)-rho(i,j,kmax))*dym1
            prgr=0.
            idpdy(i,j,kmax)=hvn(i,j,kmax)*prgr
            do k=kmax-1,kvmin_pmz(i,j),-1
               dyzu=dyzl
               dyzl=(zz(i,j+1,k)-zz(i,j,k))*dym1
               dyru=dyrl
               dyrl=(rho(i,j+1,k)-rho(i,j,k))*dym1
               dzr2=(rho(i,j+1,k+1)-rho(i,j+1,k))/(zz(i,j+1,k+1)-zz(i,j+1,k))
               dzr1=(rho(i,j  ,k+1)-rho(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
               aa=0.5*(dyrl+dyru)
               bb=0.5*(dyzl+dyzu)
               cc=0.5*(dzr2+dzr1)
               prgr=prgr+(aa-bb*cc)*0.5*(hvn(i,j,k+1)+hvn(i,j,k))
               if (k.eq.kmax-1) idpdy(i,j,kmax)=prgr*0.5*hvn(i,j,kmax)**2
               idpdy(i,j,k)=hvn(i,j,k)*prgr
            end do
         end if
      end do
   end do

   return
   end subroutine ip_blumberg_mellor_lin
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
