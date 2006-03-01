!$Id: ip_blumberg_mellor.F90,v 1.6 2006-03-01 14:45:12 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_blumberg_mellor - \label{ip-blumberg-mellor}
!
! !INTERFACE:
   subroutine ip_blumberg_mellor()
!
! !DESCRIPTION:
!
! Here, the internal part of the pressure
! gradient
! is discretised according to \cite{MELLORea94}.
! The crucial part of this term,
! which is $(\partial_x^* b)_k$ (in the case of the $u$-equation),
! is discretised
! between two vertically adjacent velocity points:
! 
! \begin{equation}\label{drhodxdiscr}
! \begin{array}{l}
! \displaystyle
! \frac12(h_{i,j,k}+h_{i,j,k+1})\left(m\,\partial_{\cal X}^*b\right)_{i,j,k} \\ \\
! \displaystyle
! \approx
! \frac12(h^u_{i,j,k}+h^u_{i,j,k+1})
! \frac{
! \frac12 (b_{i+1,j,k+1}+b_{i+1,j,k})-
! \frac12 (b_{i,j,k+1}+b_{i,j,k})}
! {\Delta x^u_{i,j}}\\ \\
! \displaystyle
! -
! \frac{z^i_{i+1,j,k}-z^i_{i,j,k}}{\Delta x^u_{i,j}}
! \left(\frac12 (b_{i+1,j,k+1}+b_{i,j,k+1})-
! \frac12 (b_{i+1,j,k}+b_{i,j,k})\right),
! \end{array}
! \end{equation}
!
! where $z^i_{i,j,k}$ is the $z$-coordinate of the interface in the T-point
! above the grid box with the index $(i,j,k)$.
!
! The discretisation of $(\partial_y^* b)_k$ for the $v$-equation is
! done accordingly.
!
! In this routine, as a first step, the interface heights are calculated
! in the T-points, in order to allow for the calculation of the 
! coordinate slopes in the U- and V-points. In a second step, the
! expression (\ref{drhodxdiscr}) equivalent formulation for the
! $y$-direction are integrated up downwards, beginning from the surface.
!
! !USES:
   use internal_pressure
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1
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
               prgr=prgr+0.5*(grdu+grdl)*0.5*(hun(i,j,k)+hun(i,j,k+1))-dxz*(rhou-rhol)
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
               prgr=prgr+0.5*(grdu+grdl)*0.5*(hvn(i,j,k)+hvn(i,j,k+1))-dyz*(rhou-rhol)
               idpdy(i,j,k)=hvn(i,j,k)*prgr
            end do
         end if
      end do
   end do

   return
   end subroutine ip_blumberg_mellor
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard, Adolf Stips and Karsten Bolding  !
!-----------------------------------------------------------------------
