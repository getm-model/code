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
   REALTYPE                  :: dxm1,dym1,zi
   REALTYPE                  :: prgr,dxzu,dxzl,dyzu,dyzl,hiu,hil
   REALTYPE                  :: dzr2,dzr1,dxru,dxrl,dyru,dyrl,aa,bb,cc
   REALTYPE,dimension(I3DFIELD) :: hw
!EOP
!-----------------------------------------------------------------------
!BOC

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
#ifndef SLICE_MODEL
   idpdy(:,:,0) = _ZERO_
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k)                                          &
!$OMP          PRIVATE(dxm1,dym1,zi)                                   &
!$OMP          PRIVATE(prgr,dxzu,dxzl,dyzu,dyzl,hiu,hil)               &
!$OMP          PRIVATE(dzr2,dzr1,dxru,dxrl,dyru,dyrl,aa,bb,cc)

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the pressure point heights are calculated in order to get the
!  interface slopes further down.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            zi = -H(i,j)
            do k=1,kmax
               hw(i,j,k-1) = _HALF_ * ( hn(i,j,k-1) + hn(i,j,k) )
               zz(i,j,k) = zi + _HALF_*hn(i,j,k)
               zi = zi + hn(i,j,k)
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
            hil=_HALF_*(hn(i,j,kmax)+hn(i+1,j,kmax))
            dxzl=(zz(i+1,j,kmax)-zz(i,j,kmax))*dxm1
            dxrl=(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            prgr=_HALF_*hil*dxrl
            idpdx(i,j,kmax)=hil*prgr
            do k=kmax-1,kumin_pmz(i,j),-1
               hiu=hil
               hil=_HALF_*(hn(i,j,k)+hn(i+1,j,k))
               dxzu=dxzl
               dxzl=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               dxru=dxrl
               dxrl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               dzr1 = ( buoy(i  ,j,k+1) - buoy(i  ,j,k) ) / hw(i  ,j,k)
               dzr2 = ( buoy(i+1,j,k+1) - buoy(i+1,j,k) ) / hw(i+1,j,k)
               aa=_HALF_*(dxrl+dxru)
               bb=_HALF_*(dxzl+dxzu)
               cc=_HALF_*(dzr1+dzr2)
               prgr=prgr+_HALF_*(hil+hiu)*(aa-bb*cc)
               idpdx(i,j,k)=hil*prgr
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1 = _ONE_/DYV
#endif
            hil=_HALF_*(hn(i,j,kmax)+hn(i,j+1,kmax))
            dyzl=(zz(i,j+1,kmax)-zz(i,j,kmax))*dym1
            dyrl=(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            prgr=_HALF_*hil*dyrl
            idpdy(i,j,kmax)=hil*prgr
            do k=kmax-1,kvmin_pmz(i,j),-1
               hiu=hil
               hil=_HALF_*(hn(i,j,k)+hn(i,j+1,k))
               dyzu=dyzl
               dyzl=(zz(i,j+1,k)-zz(i,j,k))*dym1
               dyru=dyrl
               dyrl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               dzr1 = ( buoy(i,j  ,k+1) - buoy(i  ,j,k) ) / hw(i  ,j,k)
               dzr2 = ( buoy(i,j+1,k+1) - buoy(i,j+1,k) ) / hw(i,j+1,k)
               aa=_HALF_*(dyrl+dyru)
               bb=_HALF_*(dyzl+dyzu)
               cc=_HALF_*(dzr1+dzr2)
               prgr=prgr+_HALF_*(hil+hiu)*(aa-bb*cc)
               idpdy(i,j,k)=hil*prgr
            end do
         end if
      end do
   end do
!$OMP END DO

#endif

!$OMP END PARALLEL

   return
   end subroutine ip_blumberg_mellor_lin
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
