#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ip_blumberg_mellor - \label{ip-blumberg-mellor}
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
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Adolf Stips, Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,kk,kkp
   REALTYPE                  :: dxm1,dym1,pdiff
   REALTYPE,dimension(0:1)   :: hvel,buoyvel,buoydiff
!EOP
!-----------------------------------------------------------------------
!BOC

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
#ifndef SLICE_MODEL
   idpdy(:,:,0) = _ZERO_
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k,kk,kkp)                                   &
!$OMP          PRIVATE(dxm1,dym1,pdiff)                                &
!$OMP          PRIVATE(hvel,buoyvel,buoydiff)

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the interface heights are calculated in order to get the
!  interface slopes further down.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ge. 1) then
            zz(i,j,0) = -H(i,j)
            do k=1,kmax
               zz(i,j,k) = zz(i,j,k-1) + hn(i,j,k)
            end do
         end if
      end do
   end do
!$OMP END DO


!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the u-velocity equation.

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1 = _ONE_/DXU
#endif
            kk = 0
            hvel(kk) = _HALF_ * ( hn(i,j,kmax) + hn(i+1,j,kmax) )
            buoyvel(kk) = _HALF_ * ( buoy(i,j,kmax) + buoy(i+1,j,kmax) )
            buoydiff(kk) = buoy(i+1,j,kmax) - buoy(i,j,kmax)
            pdiff =   buoyvel(kk)*( ssen(i+1,j) - ssen(i,j) ) &
                    + _HALF_*hvel(kk)*buoydiff(kk)
            idpdx(i,j,kmax) = hvel(kk)*pdiff*dxm1
            do k=kmax-1,1,-1
               kkp = kk
               kk  = 1-kk
               hvel(kk) = _HALF_ * ( hn(i,j,k) + hn(i+1,j,k) )
               buoyvel(kk) = _HALF_ * ( buoy(i,j,k) + buoy(i+1,j,k) )
               buoydiff(kk) = buoy(i+1,j,k) - buoy(i,j,k)
               pdiff =   pdiff &
                       + _QUART_*(hvel(kk)+hvel(kkp))*(buoydiff(kk)+buoydiff(kkp)) &
                       - ( zz(i+1,j,k) - zz(i,j,k) )*( buoyvel(kkp) - buoyvel(kk) )
               idpdx(i,j,k) = hvel(kk)*pdiff*dxm1
            end do
         end if
      end do
   end do
!$OMP END DO


#ifndef SLICE_MODEL

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the v-velocity equation.

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1 = _ONE_/DYV
#endif
            kk = 0
            hvel(kk) = _HALF_ * ( hn(i,j,kmax) + hn(i,j+1,kmax) )
            buoyvel(kk) = _HALF_ * ( buoy(i,j,kmax) + buoy(i,j+1,kmax) )
            buoydiff(kk) = buoy(i,j+1,kmax) - buoy(i,j,kmax)
            pdiff =   buoyvel(kk)*( ssen(i,j+1) - ssen(i,j) ) &
                    + _HALF_*hvel(kk)*buoydiff(kk)
            idpdy(i,j,kmax) = hvel(kk)*pdiff*dym1
            do k=kmax-1,1,-1
               kkp = kk
               kk  = 1-kk
               hvel(kk) = _HALF_ * ( hn(i,j,k) + hn(i,j+1,k) )
               buoyvel(kk) = _HALF_ * ( buoy(i,j,k) + buoy(i,j+1,k) )
               buoydiff(kk) = buoy(i,j+1,k) - buoy(i,j,k)
               pdiff =   pdiff &
                       + _QUART_*(hvel(kk)+hvel(kkp))*(buoydiff(kk)+buoydiff(kkp)) &
                       - ( zz(i,j+1,k) - zz(i,j,k) )*( buoyvel(kkp) - buoyvel(kk) )
               idpdy(i,j,k) = hvel(kk)*pdiff*dym1
            end do
         end if
      end do
   end do
!$OMP END DO

#endif


!$OMP END PARALLEL

   return
   end subroutine ip_blumberg_mellor
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard, Adolf Stips and Karsten Bolding  !
!-----------------------------------------------------------------------
