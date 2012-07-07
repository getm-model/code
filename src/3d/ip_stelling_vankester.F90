#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ip_stelling_vankester
!
! !INTERFACE:
   subroutine ip_stelling_vankester()
!
! !DESCRIPTION:
!
! Here, the horizontal gradients of buoyancy, $(\partial_x^* b)_k$ and
! $(\partial_y^* b)_k$, are calculated as suggested in
! Stelling and vanKester (1994).
! The horizontal gradient of buoyancy is calculated with defining
! kmax non-sloping control volumes in each water column and evaluating the
! horizontal gradients at the intersections of neighbouring control volumes.
! For each intersection, the buoyancy gradient is evaluated by linear
! interpolation of the buoyancy profile in the  neighbour column
! at the T-depth of the actual column for both directions. The minimum
! of the absolute value of the buoyancy gradient for both directions is
! used then for the internal pressure calculation.
! If both gradients point inconcistently in different directions, the
! buoyancy gradient in an intersection does not contribute to the
! internal pressure (as happens for violated hydrostatic consistency and strong
! stratification)
!
! !USES:
   use internal_pressure
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,l,kcount, rc
   REALTYPE                  :: dxm1,dym1
   REALTYPE                  :: prgr,dyz,dzz,zlm
   integer                   :: klower,kupper
   integer                   :: lnum
   REALTYPE                  :: db,dcn,dcm
   logical                   :: changed
   REALTYPE                  :: zltmp
   REALTYPE                  :: buoyplus,buoyminus
   REALTYPE                  :: zi(I3DFIELD)
   REALTYPE, POINTER         :: zx(:)
   REALTYPE, POINTER         :: zl(:)
   REALTYPE, POINTER         :: dzl(:)
   REALTYPE, POINTER         :: dzfrac(:)
   integer, POINTER          :: lvel(:)
   integer, POINTER          :: m(:)
   integer, POINTER          :: n(:)
!EOP
!-----------------------------------------------------------------------
!BOC
   lnum=2*kmax+1
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

! BJB-TODO: The zeroing of these three arrays is costly.
!  Try to reduce amount of initialization. BJB 2009-09-22.
   zz(:,:,0) = _ZERO_

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,l,kcount, rc)                                   &
!$OMP    PRIVATE(prgr,dyz,dzz,zlm,klower,kupper,lnum,db,dcn,dcm)       &
!$OMP    PRIVATE(changed,zltmp,buoyplus,buoyminus)                     &
!$OMP    PRIVATE(zx,zl,dzl,dzfrac,lvel,m,n)

! OMP-NOTE: Each thread allocates its own HEAP storage for the
!    vertical work storage:
   allocate(zx(kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (zx)'
   allocate(zl(2*kmax+1),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (zl)'
   allocate(dzl(2*kmax+1),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (dzl)'
   allocate(dzfrac(kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (dzfrac)'
   allocate(lvel(kmax+1),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (lvel)'
   allocate(m(2*kmax+1),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (m)'
   allocate(n(2*kmax+1),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (n)'

! Test de-init for funky compilers:
! NOTE: It is not necessary to initialize the 1D work arrays
!    nor the 3D array zi. BJB 2009-09-26.

! OMP-NOTE: Master thread initialize while other threads go ahead
!   and do work.
!$OMP MASTER
   idpdx(:,:,0) = _ZERO_
   idpdy(:,:,0) = _ZERO_
!$OMP END MASTER


!  First, the heights of the pressure points are calculated.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+_HALF_*hn(i,j,1)
            zi(i,j,1)=-H(i,j)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+_HALF_*(hn(i,j,k-1)+hn(i,j,k))
               zi(i,j,k)=zi(i,j,k-1)+hn(i,j,k-1)
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
         n(:)=0
         m(:)=0
         zl(:)=_ZERO_

         if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dxm1=_ONE_/DXU
#endif
            ! merge zl vector from zi's
            zx(1)=-HU(i,j)+_HALF_*hun(i,j,1) ! zx defined on u-points
            zl(1)=zi(i,j,1)
            zl(2)=zi(i+1,j,1)
            do k=2,kmax
               zx(k)=zx(k-1)+_HALF_*(hun(i,j,k-1)+hun(i,j,k))
               zl(k*2-1)=zi(i,j,k)
               zl(k*2)=zi(i+1,j,k)
            end do
            zl(lnum)=ssun(i,j)

            ! sort the zl vector
            do
               changed=.false.
               do l=1,lnum-1
                  if (zl(l) .gt. zl(l+1)) then
                      zltmp=zl(l+1)
                      zl(l+1)=zl(l)
                      zl(l)=zltmp
                      changed=.true.
                  end if
               end do
               if (.not.changed) exit
            end do

            ! calc dzl, find m,n and idpdx-points
            m(lnum)=kmax !left cell number
            n(lnum)=kmax !right cell number
            kcount=kmax
            lvel(kmax+1)=lnum
            do l=lnum-1,1,-1
               dzl(l)=zl(l+1)-zl(l)
               m(l)=m(l+1)
               n(l)=n(l+1)
               if (zl(l+1) .le. zi(i  ,j,m(l))) m(l)=max(1,m(l)-1)
               if (zl(l+1) .le. zi(i+1,j,n(l))) n(l)=max(1,n(l)-1)

               if (kcount .gt. 0) then
                  if (zx(kcount) .gt. zl(l)) then
                     lvel(kcount)=l
                     dzfrac(kcount)=(zx(kcount)-zl(l))
                     kcount = kcount-1
                  end if
               end if
            end do

            ! integrate internal pressure
            prgr=_ZERO_
            db=_ZERO_
            do k=kmax,1,-1
               do l=lvel(k+1)-1,lvel(k),-1
                  ! calc b-gradient
                  zlm=_HALF_*dzl(l)+zl(l)

                  if (zz(i+1,j,n(l)).ge.zz(i,j,kmax)) then
                           buoyminus=buoy(i,j,kmax)
                  else
                    if (zz(i+1,j,n(l)).le.zz(i,j,1)) then
                            buoyminus=buoy(i,j,1)
                    else
                      if (zz(i+1,j,n(l)).ge.zz(i,j,m(l))) then
                          klower=m(l)
                          kupper=m(l)+1
                      else
                          klower=m(l)-1
                          kupper=m(l)
                      end if
                      dzz=zz(i,j,kupper)-zz(i,j,klower)
                      buoyminus=(zz(i+1,j,n(l))-zz(i,j,klower))/dzz * &
                            buoy(i,j,kupper) + &
                            (zz(i,j,kupper)-zz(i+1,j,n(l)))/dzz * &
                            buoy(i,j,klower)
                    end if
                  end if

                  if (zz(i,j,m(l)).ge.zz(i+1,j,kmax)) then
                           buoyplus=buoy(i+1,j,kmax)
                  else
                    if (zz(i,j,m(l)).le.zz(i+1,j,1)) then
                           buoyplus=buoy(i+1,j,1)
                    else
                      if (zz(i,j,m(l)).ge.zz(i+1,j,n(l))) then
                          klower=n(l)
                          kupper=n(l)+1
                      else
                          klower=n(l)-1
                          kupper=n(l)
                      end if
                      dzz=zz(i+1,j,kupper)-zz(i+1,j,klower)
                      buoyplus=(zz(i,j,m(l))-zz(i+1,j,klower))/dzz *  &
                           buoy(i+1,j,kupper) + &
                           (zz(i+1,j,kupper)-zz(i,j,m(l)))/dzz *  &
                           buoy(i+1,j,klower)
                    end if
                  end if

                  dcm=buoy(i+1,j,n(l))-buoyminus
                  dcn=buoyplus-buoy(i,j,m(l))

                  if (dcn*dcm .lt. _ZERO_) then
                     db=_ZERO_
                  else
                     if (dcn.ge._ZERO_) then
                         db=min(dcn,dcm)
                     else
                         db=max(dcn,dcm)
                     end if
                  end if

                  !calc p-gradient (even in the dzl of u-point)
                  prgr = prgr + dzl(l)*db*dxm1
               end do
               idpdx(i,j,k) = hun(i,j,k)*(prgr - &
                                dzfrac(k)*db*dxm1)
            end do
         end if
      end do
   end do
!$OMP END DO

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         n(:)=0
         m(:)=0
         zl(:)=_ZERO_

         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dym1 = _ONE_/DYV
#endif
            ! merge zl vector from zi's
            zx(1)=-HV(i,j)+_HALF_*hvn(i,j,1) ! zx defined on u-points
            zl(1)=zi(i,j,1)
            zl(2)=zi(i,j+1,1)
            do k=2,kmax
               zx(k)=zx(k-1)+_HALF_*(hvn(i,j,k-1)+hvn(i,j,k))
               zl(k*2-1)=zi(i,j,k)
               zl(k*2)=zi(i,j+1,k)
            end do
            zl(lnum)=ssvn(i,j)

            ! sort the zl vector
            do
               changed=.false.
               do l=1,lnum-1
                  if (zl(l) .gt. zl(l+1)) then
                      zltmp=zl(l+1)
                      zl(l+1)=zl(l)
                      zl(l)=zltmp
                      changed=.true.
                  end if
               end do
               if (.not.changed) exit
            end do

            ! calc dzl, find m,n and idpdx-points
            m(lnum)=kmax !left cell number
            n(lnum)=kmax !right cell number
            kcount=kmax
            lvel(kmax+1)=lnum
            do l=lnum-1,1,-1
               dzl(l)=zl(l+1)-zl(l)
               m(l)=m(l+1)
               n(l)=n(l+1)
               if (zl(l+1) .le. zi(i  ,j,m(l))) m(l)=max(1,m(l)-1)
               if (zl(l+1) .le. zi(i,j+1,n(l))) n(l)=max(1,n(l)-1)

               if (kcount .gt. 0) then
                  if (zx(kcount) .gt. zl(l)) then
                     lvel(kcount)=l
                     dzfrac(kcount)=(zx(kcount)-zl(l))
                     kcount = kcount-1
                  end if
               end if
            end do

            ! integrate internal pressure
            prgr=_ZERO_
            do k=kmax,1,-1
               do l=lvel(k+1)-1,lvel(k),-1
                  ! calc b-gradient
                  zlm=_HALF_*dzl(l)+zl(l)

                  if (zz(i,j+1,n(l)).ge.zz(i,j,kmax)) then
                           buoyminus=buoy(i,j,kmax)
                  else
                    if (zz(i,j+1,n(l)).le.zz(i,j,1)) then
                            buoyminus=buoy(i,j,1)
                    else
                      if (zz(i,j+1,n(l)).ge.zz(i,j,m(l))) then
                          klower=m(l)
                          kupper=m(l)+1
                      else
                          klower=m(l)-1
                          kupper=m(l)
                      end if
                      dzz=zz(i,j,kupper)-zz(i,j,klower)
                      buoyminus=(zz(i,j+1,n(l))-zz(i,j,klower))/dzz * &
                            buoy(i,j,kupper) + &
                            (zz(i,j,kupper)-zz(i,j+1,n(l)))/dzz * &
                            buoy(i,j,klower)
                    end if
                  end if

                  if (zz(i,j,m(l)).ge.zz(i,j+1,kmax)) then
                           buoyplus=buoy(i,j+1,kmax)
                  else
                    if (zz(i,j,m(l)).le.zz(i,j+1,1)) then
                           buoyplus=buoy(i,j+1,1)
                    else
                      if (zz(i,j,m(l)).ge.zz(i,j+1,n(l))) then
                          klower=n(l)
                          kupper=n(l)+1
                      else
                          klower=n(l)-1
                          kupper=n(l)
                      end if
                      dzz=zz(i,j+1,kupper)-zz(i,j+1,klower)
                      buoyplus=(zz(i,j,m(l))-zz(i,j+1,klower))/dzz *  &
                           buoy(i,j+1,kupper) + &
                           (zz(i,j+1,kupper)-zz(i,j,m(l)))/dzz *  &
                           buoy(i,j+1,klower)
                    end if
                  end if

                  dcm=buoy(i,j+1,n(l))-buoyminus
                  dcn=buoyplus-buoy(i,j,m(l))

                  if (dcn*dcm .lt. _ZERO_) then
                     db=_ZERO_
                  else
                     if (dcn.ge._ZERO_) db=min(dcn,dcm)
                     if (dcn.le._ZERO_) db=max(dcn,dcm)
                  end if

                  !calc p-gradient (even in the dzl of u-point)
                  prgr = prgr + dzl(l)*db*dym1
               end do
               idpdy(i,j,k) = hvn(i,j,k)*(prgr - &
                                dzfrac(k)*db*dym1)
            end do
         end if
      end do
   end do
!$OMP END DO

! Each thread must deallocate its own HEAP storage:
   deallocate(zx,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (zx)'
   deallocate(zl,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (zl)'
   deallocate(dzl,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (dzl)'
   deallocate(dzfrac,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (dzfrac)'
   deallocate(lvel,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (lvel)'
   deallocate(m,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (m)'
   deallocate(n,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (n)'

!$OMP END PARALLEL
   return
   end subroutine ip_stelling_vankester
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
