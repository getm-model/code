!$Id: ip_blumberg_mellor_lin.F90,v 1.5 2007-02-23 12:20:36 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_blumberg_mellor_lin
!
! !INTERFACE:
   subroutine ip_blumberg_mellor_lin()
!
! !DESCRIPTION:
!
! Here, the internal pressure gradient calculation is carried out on the
! basis of the same buoyancy stencil than in the method according
! to \cite{MELLORea94} (see routine {\tt ip\_blumberg\_mellor}), 
! but in such a way that the pressure gradient numerically vanishes for
! linear stratification without horizontal gradients.
!
! \begin{equation}\label{drhodxdiscr_lin}
! \begin{array}{l}
! \displaystyle
! \frac12(h_{i,j,k}+h_{i,j,k+1})\left(m\,\partial_{\cal X}^*b\right)_{i,j,k} \\
! \\
! \displaystyle
! \approx
! \frac12(h^u_{i,j,k}+h^u_{i,j,k+1})
! \Bigg[\frac{
! \frac12 (b_{i+1,j,k+1}+b_{i+1,j,k})-
! \frac12 (b_{i,j,k+1}+b_{i,j,k})}
! {\Delta x^u_{i,j}}\\ \\
! \displaystyle
! \qquad -
! \frac12\left(\frac{
! \frac12 (z^c_{i+1,j,k+1}+z^c_{i+1,j,k})-
! \frac12 (z^c_{i,j,k+1}+z^c_{i,j,k})}
! {\Delta x^u_{i,j}}\right) \\ \\
! \displaystyle
! \qquad\qquad 
! \frac12\left(
! \frac{b_{i+1,j,k+1}-b_{i+1,j,k}}{z^c_{i+1,j,k+1}-z^c_{i+1,j,k}}+
! \frac{b_{i,j,k+1}-b_{i,j,k}}{z^c_{i,j,k+1}-z^c_{i,j,k}}
! \right)\Bigg],
! \end{array}
! \end{equation}
!
! where $z^c_{i,j,k}$ is the $z$-coordinate of the centre of
! the grid box with the index $(i,j,k)$.
!
! The discretisation of $(\partial_y^* b)_k$ for the $v$-equation is
! done accordingly.
!
! !USES:
   use internal_pressure
   use variables_3d, only: kumin_pmz,kvmin_pmz
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
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
            dxrl=(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            prgr=0.
            do k=kmax-1,kumin_pmz(i,j),-1
               dxzu=dxzl
               dxzl=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               dxru=dxrl
               dxrl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               dzr2=(buoy(i+1,j,k+1)-buoy(i+1,j,k))/(zz(i+1,j,k+1)-zz(i+1,j,k))
               dzr1=(buoy(i  ,j,k+1)-buoy(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
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
           dyrl=(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            prgr=0.
            idpdy(i,j,kmax)=hvn(i,j,kmax)*prgr
            do k=kmax-1,kvmin_pmz(i,j),-1
               dyzu=dyzl
               dyzl=(zz(i,j+1,k)-zz(i,j,k))*dym1
               dyru=dyrl
               dyrl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               dzr2=(buoy(i,j+1,k+1)-buoy(i,j+1,k))/(zz(i,j+1,k+1)-zz(i,j+1,k))
               dzr1=(buoy(i,j  ,k+1)-buoy(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
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
