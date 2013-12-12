#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ip_song_wright
!
! !INTERFACE:
   subroutine ip_song_wright()
!
! !DESCRIPTION:
!
! Here, the pressure gradient is calculating according to an energy-conserving
! method suggested by \cite{SONG98}, which for the pressure gradient in
! $x$-direction looks as:
!
! \begin{equation}\label{drhodxdiscrSONG}
! \begin{array}{l}
! \displaystyle
! \frac12(h_{i,j,k}+h_{i,j,k+1})\left(m\,\partial_{\cal X}^*b\right)_{i,j,k} \\ \\
! \displaystyle
! \approx
! \frac{
! \frac14 (b_{i+1,j,k+1}+b_{i+1,j,k})(h^c_{i+1,j,k+1}+h^c_{i+1,j,k})-
! \frac14 (b_{i,j,k+1}+b_{i,j,k})(h^c_{i,j,k+1}+h^c_{i,j,k})}
! {\Delta x^u_{i,j}}\\ \\
! \displaystyle
! \qquad-
! \Bigg[\frac12 (b_{i+1,j,k+1}+b_{i,j,k+1})
! \frac{z^c_{i+1,j,k+1}-z^c_{i,j,k+1}}{\Delta x^u_{i,j}}\\ \\
! \displaystyle
! \qquad\qquad -
! \frac12 (b_{i+1,j,k}+b_{i,j,k})
! \frac{z^c_{i+1,j,k}-z^c_{i,j,k}}{\Delta x^u_{i,j}}\Bigg],
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
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1,zi
   REALTYPE                  :: grdl,grdu,buoyl,buoyu,prgr,dxz,dyz,hi
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
!$OMP          PRIVATE(grdl,grdu,buoyl,buoyu,prgr,dxz,dyz,hi)

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the interface heights are calculated in order to get the
!  interface slopes further down.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            zi = -H(i,j)
            do k=1,kmax
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
            hi=_HALF_*(hn(i,j,kmax)+hn(i+1,j,kmax))
            grdl=_HALF_*hi*(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            buoyl=_HALF_*(buoy(i,j,kmax)+buoy(i+1,j,kmax))       &
                     *(zz(i+1,j,kmax)- zz(i,j,kmax))*dxm1
            prgr=grdl
            idpdx(i,j,kmax)=hi*prgr
            do k=kmax-1,1,-1
!              KK-TODO: for what do we need grdu?
               grdu=grdl
               grdl=(_HALF_*(buoy(i+1,j,k)+buoy(i+1,j,k+1))               &
                        *_HALF_*(hn(i+1,j,k)+hn(i+1,j,k+1))             &
                    -_HALF_*(buoy(i,j,k)+buoy(i,j,k+1))                   &
                        *_HALF_*(hn(i,j,k)+hn(i,j,k+1)) )*dxm1
               buoyu=buoyl
               buoyl=_HALF_*(buoy(i,j,k)+buoy(i+1,j,k))*(zz(i+1,j,k)-zz(i,j,k))*dxm1
               prgr=prgr+grdl-(buoyu-buoyl)
               idpdx(i,j,k)=_HALF_*(hn(i,j,k)+hn(i+1,j,k))*prgr
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the v-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dym1 = _ONE_/DYV
#endif
            hi=_HALF_*(hn(i,j,kmax)+hn(i,j+1,kmax))
            grdl=_HALF_*hi*(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            buoyl=_HALF_*(buoy(i,j,kmax)+buoy(i,j+1,kmax))     &
                     *(zz(i,j+1,kmax)- zz(i,j,kmax))*dxm1
            prgr=grdl
            idpdy(i,j,kmax)=hi*prgr
            do k=kmax-1,1,-1
!              KK-TODO: for what do we need grdu?
               grdu=grdl
               grdl=(_HALF_*(buoy(i,j+1,k)+buoy(i,j+1,k+1))               &
                        *_HALF_*(hn(i,j+1,k)+hn(i,j+1,k+1))             &
                    -_HALF_*(buoy(i,j,k)+buoy(i,j,k+1))                   &
                        *_HALF_*(hn(i,j,k)+hn(i,j,k+1)) )*dym1
               buoyu=buoyl
               buoyl=_HALF_*(buoy(i,j,k)+buoy(i,j+1,k))*(zz(i,j+1,k)-zz(i,j,k))*dym1
               prgr=prgr+grdl-(buoyu-buoyl)
               idpdy(i,j,k)=_HALF_*(hn(i,j,k)+hn(i,j+1,k))*prgr
            end do
         end if
      end do
   end do
!$OMP END DO

#endif

!$OMP END PARALLEL

   return
   end subroutine ip_song_wright
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
