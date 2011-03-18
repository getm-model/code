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
!$ use omp_lib
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

! BJB-TODO: The zeroing of these three arrays is costly.
!  Try to reduce amount of initialization. BJB 2009-09-22.
   zz(:,:,0) = _ZERO_

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,prgr,dxzu,dxzl,dyzu,dyzl)                       &
!$OMP    PRIVATE(dzr2,dzr1,dxru,dxrl,dyru,dyrl,aa,bb,cc)

!$OMP MASTER
   idpdx(:,:,0) = _ZERO_
   idpdy(:,:,0) = _ZERO_
!$OMP END MASTER

!  First, the pressure point heights are calculated in order to get the
!  interface slopes further down.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+_HALF_*hn(i,j,1)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+_HALF_*(hn(i,j,k-1)+hn(i,j,k))
            end do
         end if
      end do
   end do
!$OMP END DO

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the u-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            dxzl=(zz(i+1,j,kmax)-zz(i,j,kmax))*dxm1
            dxrl=(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            prgr=dxrl*_HALF_*hun(i,j,kmax)
            idpdx(i,j,kmax)=hun(i,j,kmax)*prgr
            do k=kmax-1,kumin_pmz(i,j),-1
               dxzu=dxzl
               dxzl=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               dxru=dxrl
               dxrl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               dzr2=(buoy(i+1,j,k+1)-buoy(i+1,j,k))/(zz(i+1,j,k+1)-zz(i+1,j,k))
               dzr1=(buoy(i  ,j,k+1)-buoy(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
               aa=_HALF_*(dxrl+dxru)
               bb=_HALF_*(dxzl+dxzu)
               cc=_HALF_*(dzr2+dzr1)
               prgr=prgr+(aa-bb*cc)*_HALF_*(hun(i,j,k+1)+hun(i,j,k))
               idpdx(i,j,k)=hun(i,j,k)*prgr
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1 = _ONE_/DYV
#endif
            dyzl=(zz(i,j+1,kmax)-zz(i,j,kmax))*dym1
            dyrl=(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            prgr=dyrl*_HALF_*hun(i,j,kmax)
            idpdy(i,j,kmax)=hvn(i,j,kmax)*prgr
            do k=kmax-1,kvmin_pmz(i,j),-1
               dyzu=dyzl
               dyzl=(zz(i,j+1,k)-zz(i,j,k))*dym1
               dyru=dyrl
               dyrl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               dzr2=(buoy(i,j+1,k+1)-buoy(i,j+1,k))/(zz(i,j+1,k+1)-zz(i,j+1,k))
               dzr1=(buoy(i,j  ,k+1)-buoy(i  ,j,k))/(zz(i  ,j,k+1)-zz(i  ,j,k))
               aa=_HALF_*(dyrl+dyru)
               bb=_HALF_*(dyzl+dyzu)
               cc=_HALF_*(dzr2+dzr1)
               prgr=prgr+(aa-bb*cc)*_HALF_*(hvn(i,j,k+1)+hvn(i,j,k))
               idpdy(i,j,k)=hvn(i,j,k)*prgr
            end do
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

   return
   end subroutine ip_blumberg_mellor_lin
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
