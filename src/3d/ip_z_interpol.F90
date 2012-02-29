#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_z_interpol
!
! !INTERFACE:
   subroutine ip_z_interpol()
!
! !DESCRIPTION:
!
! Here, the horizontal gradients of buoyancy, $(\partial_x^* b)_k$ and
! $(\partial_y^* b)_k$, are directly calculated in $z$-coordinates by
! linearly interpolating the buoyancies in the vertical to the
! evaluation point (which is the interface vertically located between
! the velocity points). In the case that extrapolations become
! necessary near the sloping surface (or more likely) near the sloping
! bottom, then the last regular buoyancy value (surface value or bottom
! value) is used.
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
   integer                   :: i,j,k,kplus,kminus, rc
   REALTYPE                  :: dxm1,dym1,zi
   REALTYPE                  :: grdl,grdu,buoyl,prgr,dxz,dyz
   REALTYPE,dimension(:),POINTER :: zx,hi
   REALTYPE                  :: buoyplus,buoyminus
!EOP
!-----------------------------------------------------------------------
!BOC

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
#ifndef SLICE_MODEL
   idpdy(:,:,0) = _ZERO_
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k,kplus,kminus, rc)                         &
!$OMP          PRIVATE(dxm1,dym1,zi)                                   &
!$OMP          PRIVATE(grdl,grdu,buoyl,prgr,dxz,dyz)                   &
!$OMP          PRIVATE(zx,hi)                                          &
!$OMP          PRIVATE(buoyplus,buoyminus)

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

! OMP-NOTE: Each thread allocates its own HEAP storage for the
!    vertical work storage:
   allocate(zx(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_z_interpol: Error allocating memory (zx)'
   allocate(hi(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_z_interpol: Error allocating memory (hi)'
   hi(0) = _ZERO_

!  First, the heights of the pressure points are calculated.
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
            zx(0) = -HU(i,j) ! zx defined on u-points
            do k=1,kmax
               hi(k) = _HALF_ * ( hn(i,j,k) + hn(i+1,j,k) )
               zx(k) = zx(k-1) + _HALF_*(hi(k-1)+hi(k))
            end do
            grdl=_HALF_*hi(kmax)*(buoy(i+1,j,kmax)-buoy(i,j,kmax))*dxm1
            buoyl=_HALF_*(buoy(i,j,kmax)+buoy(i+1,j,kmax))
            prgr=grdl
            idpdx(i,j,kmax)=hi(kmax)*prgr
            do k=kmax-1,1,-1
               grdu=grdl
               do kplus=kmax,1,-1  ! Find neighboring index to east
                  if (zz(i+1,j,kplus) .le. zx(k)) EXIT
               end do
               do kminus=kmax,1,-1 ! Find neighboring index to west
                  if (zz(i,j,kminus) .le. zx(k)) EXIT
               end do
               if (kplus .eq. kmax) buoyplus=buoy(i+1,j,kplus)
               if (kplus .lt. kmax .and. kplus .gt. 1) then ! interpolate
                  buoyplus=((zx(k)-zz(i+1,j,kplus))*buoy(i+1,j,kplus+1)+ &
                          (zz(i+1,j,kplus+1)-zx(k))*buoy(i+1,j,kplus))/ &
                          (_HALF_*(hn(i+1,j,kplus+1)+hn(i+1,j,kplus)))
               end if
               if (kminus .eq. kmax) buoyminus=buoy(i,j,kminus)
               if ((kminus .lt. kmax) .and. (kminus .gt. 1)) then ! interpolate
                  buoyminus=((zx(k)-zz(i,j,kminus))*buoy(i,j,kminus+1)+ &
                          (zz(i,j,kminus+1)-zx(k))*buoy(i,j,kminus))/  &
               (_HALF_*(hn(i,j,kminus+1)+hn(i,j,kminus)))
               end if
               if (zx(k) .gt. max(-H(i+1,j),-H(i,j))) then
                  grdl=_HALF_*hi(k)*(buoyplus-buoyminus)*dxm1
               else
                  grdl= _ZERO_
               end if
               prgr=prgr+grdu+grdl
               idpdx(i,j,k)=hi(k)*prgr
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
            zx(0) = -HV(i,j) ! zx defined on v-points
            do k=1,kmax
               hi(k) = _HALF_ * ( hn(i,j,k) + hn(i,j+1,k) )
               zx(k) = zx(k-1) + _HALF_*(hi(k-1)+hi(k))
            end do
            grdl=_HALF_*hi(kmax)*(buoy(i,j+1,kmax)-buoy(i,j,kmax))*dym1
            buoyl=_HALF_*(buoy(i,j+1,kmax)+buoy(i,j,kmax))
            prgr=grdl
            idpdy(i,j,kmax)=hi(kmax)*prgr
            do k=kmax-1,1,-1
               grdu=grdl
               do kplus=kmax,1,-1  ! Find neighboring index to north
                  if (zz(i,j+1,kplus) .le. zx(k)) EXIT
               end do
               do kminus=kmax,1,-1 ! Find neighboring index to south
                  if (zz(i,j,kminus) .le. zx(k)) EXIT
               end do
               if (kplus .eq. kmax) buoyplus=buoy(i,j+1,kplus)
               if ((kplus .lt. kmax) .and. (kplus .gt. 1)) then
                  buoyplus=((zx(k)-zz(i,j+1,kplus))*buoy(i,j+1,kplus+1)+ &
                          (zz(i,j+1,kplus+1)-zx(k))*buoy(i,j+1,kplus))/ &
                          (_HALF_*(hn(i,j+1,kplus+1)+hn(i,j+1,kplus)))
               end if
               if (kminus .eq. kmax) buoyminus=buoy(i,j,kminus)
               if ((kminus .lt. kmax) .and. (kminus .gt. 1)) then
                  buoyminus=((zx(k)-zz(i,j,kminus))*buoy(i,j,kminus+1)+  &
                           (zz(i,j,kminus+1)-zx(k))*buoy(i,j,kminus))/  &
                           (_HALF_*(hn(i,j,kminus+1)+hn(i,j,kminus)))
               end if
               if (zx(k).gt.max(-H(i,j+1),-H(i,j))) then
                  grdl=_HALF_*hi(k)*(buoyplus-buoyminus)*dym1
               else
                  grdl= _ZERO_
               end if
               prgr=prgr+grdu+grdl
               idpdy(i,j,k)=hi(k)*prgr
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

#endif

! Each thread must deallocate its own HEAP storage:
   deallocate(zx,stat=rc)
   if (rc /= 0) stop 'ip_z_interpol: Error deallocating memory (zx)'
   deallocate(hi,stat=rc)
   if (rc /= 0) stop 'ip_z_interpol: Error deallocating memory (hi)'

!$OMP END PARALLEL

   return
   end subroutine ip_z_interpol
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
