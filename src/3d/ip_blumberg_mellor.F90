!$Id: ip_blumberg_mellor.F90,v 1.10 2007-07-12 10:26:00 kbk Exp $
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
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Adolf Stips, Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1
   REALTYPE                  :: grdl,grdu,buoyl,buoyu,prgr,dxz,dyz
!EOP
!-----------------------------------------------------------------------
!BOC
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the interface heights are calculated in order to get the
!  interface slopes further down.
   do j=jmin,jmax+1
      do i=imin,imax+1
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
   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            grdl=(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            buoyl=0.5*(buoy(i+1,j,kmax)+buoy(i,j,kmax))
            prgr=grdl*0.5*hun(i,j,kmax)
            idpdx(i,j,kmax)=hun(i,j,kmax)*prgr
            do k=kmax-1,1,-1
               grdu=grdl
               grdl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               buoyu=buoyl
               buoyl=0.5*(buoy(i+1,j,k)+buoy(i,j,k))
               dxz=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               prgr=prgr+0.5*(grdu+grdl)*0.5*(hun(i,j,k)+hun(i,j,k+1))-dxz*(buoyu-buoyl)
               idpdx(i,j,k)=hun(i,j,k)*prgr
            end do
         end if
      end do
   end do

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the v-velocity equation.
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1 = _ONE_/DYV
#endif
            grdl=(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            buoyl=0.5*(buoy(i,j+1,kmax)+buoy(i,j,kmax))
            prgr=grdl*0.5*hvn(i,j,kmax)
            idpdy(i,j,kmax)=hvn(i,j,kmax)*prgr
            do k=kmax-1,1,-1
               grdu=grdl
               grdl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               buoyu=buoyl
               buoyl=0.5*(buoy(i,j+1,k)+buoy(i,j,k))
               dyz=(zz(i,j+1,k)-zz(i,j,k))*dym1
               prgr=prgr+0.5*(grdu+grdl)*0.5*(hvn(i,j,k)+hvn(i,j,k+1))-dyz*(buoyu-buoyl)
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
