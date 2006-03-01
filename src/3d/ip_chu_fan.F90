!$Id: ip_chu_fan.F90,v 1.3 2006-03-01 14:45:12 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_chu_fan 
!
! !INTERFACE:
   subroutine ip_chu_fan()
!
! !DESCRIPTION:
!   
! This routine calculates the internal pressure gradient based on the
! classical approach by \cite{MELLORea94}, extended by the
! hydrostatic extension by \cite{CHUea03}.
!
! !USES:
   use internal_pressure
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1,x,y,x1,y1,hc
   REALTYPE                  :: grdl,grdu,rhol,rhou,prgr,dxz,dyz
!EOP
!-----------------------------------------------------------------------
!BOC
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the interface heights are calculated in order to get the
!  interface slopes further down.
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

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the u-velocity equation.
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
               hc=(zz(i+1,j,k)-zz(i+1,j,k+1))**2*(rho(i+1,j,k)-rho(i+1,j,k+1))
               hc=hc-(zz(i,j,k)-zz(i,j,k+1))**2*(rho(i,j,k)-rho(i,j,k+1)) 
               hc=hc+(zz(i+1,j,k+1)-zz(i,j,k+1))**2*(rho(i+1,j,k+1)-rho(i,j,k+1)) 
               hc=hc-(zz(i+1,j,k)-zz(i,j,k))**2*(rho(i+1,j,k)-rho(i,j,k)) 
               hc=hc*dxm1*0.1666666667/(zz(i,j,k)+zz(i+1,j,k)-zz(i,j,k+1)-zz(i+1,j,k+1))
               prgr=prgr+(grdu+grdl+hc)*0.5*(hun(i,j,k)+hun(i,j,k+1))-dxz*(rhou-rhol)
               idpdx(i,j,k)=hun(i,j,k)*prgr
            end do
         end if
      end do
   end do

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the v-velocity equation.
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
               hc=(zz(i,j+1,k)-zz(i,j+1,k+1))**2*(rho(i,j+1,k)-rho(i,j+1,k+1))
               hc=hc-(zz(i,j,k)-zz(i,j,k+1))**2*(rho(i,j,k)-rho(i,j,k+1)) 
               hc=hc+(zz(i,j+1,k+1)-zz(i,j,k+1))**2*(rho(i,j+1,k+1)-rho(i,j,k+1)) 
               hc=hc-(zz(i,j+1,k)-zz(i,j,k))**2*(rho(i,j+1,k)-rho(i,j,k)) 
               hc=hc*dym1*0.1666666667/(zz(i,j,k)+zz(i,j+1,k)-zz(i,j,k+1)-zz(i,j+1,k+1))
               prgr=prgr+(grdu+grdl+hc)*0.5*(hvn(i,j,k)+hvn(i,j,k+1))-dyz*(rhou-rhol)
               idpdy(i,j,k)=hvn(i,j,k)*prgr
            end do
         end if
      end do
   end do

   return
   end subroutine ip_chu_fan
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
