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
! Note (KK): The hydrostatic correction approach of \cite{CHUea03} is
! based on pressure-jacobian scheme and probably cannot be mixed with
! density-jacobian scheme of \cite{MELLORea94}.
!
! !USES:
   use internal_pressure
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding & Adolf Stips
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dxm1,dym1,x,y,x1,y1,hc
   REALTYPE                  :: grdl,grdu,buoyl,buoyu,prgr,dxz,dyz,hil,hiu
   REALTYPE, PARAMETER       :: SIXTH=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC

! OMP-NOTE: The initialization and OMP implementation in this routine
!   is not tested, as the initialization states that the present
!   method is "Not working, use other internal pressure gradient scheme"
!   BJB 2009-09-24.

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
   idpdy(:,:,0) = _ZERO_

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k)                                          &
!$OMP          PRIVATE(dxm1,dym1,x,y,x1,y1,hc)                         &
!$OMP          PRIVATE(grdl,grdu,buoyl,buoyu,prgr,dxz,dyz,hil,hiu)

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
   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            hil=_HALF_*(hn(i,j,kmax)+hn(i+1,j,kmax))
            grdl=(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            buoyl=_HALF_*(buoy(i,j,kmax)+buoy(i+1,j,kmax))
            prgr=grdl
!           here prgr has units of buoy/m, idpdx has correct units
            idpdx(i,j,kmax)=hil*prgr*_HALF_*hil
            do k=kmax-1,1,-1
               hiu=hil
               hil=_HALF_*(hn(i,j,k)+hn(i+1,j,k))
               grdu=grdl
               grdl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               buoyu=buoyl
               buoyl=_HALF_*(buoy(i,j,k)+buoy(i+1,j,k))
               dxz=(zz(i+1,j,k)-zz(i,j,k))*dxm1
               hc =    - hn(i+1,j,k+1)**2 * ( buoy(i+1,j,k+1) - buoy(i+1,j,k) )
               hc = hc + hn(i  ,j,k+1)**2 * ( buoy(i  ,j,k+1) - buoy(i  ,j,k) )
               hc=hc+(zz(i+1,j,k+1)-zz(i,j,k+1))**2*(buoy(i+1,j,k+1)-buoy(i,j,k+1))
               hc=hc-(zz(i+1,j,k)-zz(i,j,k))**2*(buoy(i+1,j,k)-buoy(i,j,k))
               hc=-hc*dxm1*SIXTH/(hn(i,j,k+1)+hn(i+1,j,k+1))
!              hc has units of buoy, grd[u|l] has units of buoy/m, mixture of units for prgr
               prgr=prgr+(grdl+grdu+hc)*_HALF_*(hil+hiu)-dxz*(buoyu-buoyl)
               idpdx(i,j,k)=hil*prgr
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

!  Calculation of layer integrated internal pressure gradient as it
!  appears on the right hand side of the v-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1 = _ONE_/DYV
#endif
            hil=_HALF_*(hn(i,j,kmax)+hn(i,j+1,kmax))
            grdl=(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            buoyl=_HALF_*(buoy(i,j,kmax)+buoy(i,j+1,kmax))
            prgr=grdl
            idpdy(i,j,kmax)=hil*prgr*_HALF_*hil
            do k=kmax-1,1,-1
               hiu=hil
               hil=_HALF_*(hn(i,j,k)+hn(i,j+1,k))
               grdu=grdl
               grdl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               buoyu=buoyl
               buoyl=_HALF_*(buoy(i,j,k)+buoy(i,j+1,k))
               dyz=(zz(i,j+1,k)-zz(i,j,k))*dym1
               hc =    - hn(i,j+1,k+1)**2 * ( buoy(i,j+1,k+1) - buoy(i,j+1,k) )
               hc = hc + hn(i,j  ,k+1)**2 * ( buoy(i,j  ,k+1) - buoy(i,j  ,k) )
               hc=hc+(zz(i,j+1,k+1)-zz(i,j,k+1))**2*(buoy(i,j+1,k+1)-buoy(i,j,k+1))
               hc=hc-(zz(i,j+1,k)-zz(i,j,k))**2*(buoy(i,j+1,k)-buoy(i,j,k))
               hc=-hc*dym1*SIXTH/(hn(i,j,k+1)+hn(i,j+1,k+1))
               prgr=prgr+(grdl+grdu+hc)*_HALF_*(hil+hiu)-dyz*(buoyu-buoyl)
               idpdy(i,j,k)=hil*prgr
            end do
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL

   return
   end subroutine ip_chu_fan
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
