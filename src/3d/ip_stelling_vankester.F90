#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_stelling_vankester
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
   integer                           :: i,j,k,kl,kr,l,lmin,lmax,lextr,rc
   REALTYPE                          :: dxm1,dym1
   REALTYPE                          :: pdiff,buoydiff
   REALTYPE                          :: zvel,zm,zn,dcm,dcn,corr,weight
   REALTYPE,dimension(I3DFIELD)      :: zi
   REALTYPE,dimension(:),allocatable :: zl
   integer,dimension(:),allocatable  :: m,n
!EOP
!-----------------------------------------------------------------------
!BOC

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
   idpdy(:,:,0) = _ZERO_

   lextr = 2*kmax+1

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k,kl,kr,l,lmin,lmax,rc)                     &
!$OMP          PRIVATE(dxm1,dym1)                                      &
!$OMP          PRIVATE(pdiff,buoydiff)                                 &
!$OMP          PRIVATE(zvel,zm,zn,dcm,dcn,corr,weight)                 &
!$OMP          PRIVATE(zl,m,n)

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

! OMP-NOTE: Each thread allocates its own HEAP storage for the
!    vertical work storage:
   allocate(zl(0:lextr),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (zl)'
   allocate(m(0:lextr),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (m)'
   allocate(n(0:lextr),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_stelling_vankester: Error allocating memory (n)'

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ge. 1) then
            zi(i,j,0) = -H(i,j)
            do k=1,kmax
               zz(i,j,k) = zi(i,j,k-1) + _HALF_*hn(i,j,k)
               zi(i,j,k) = zi(i,j,k-1) + hn(i,j,k)
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
            dxm1=_ONE_/DXU
#endif
            l = 0
            kl = 0
            kr = 0
            lmax = -1 ! should cause segmentation fault
            do while(kl.le.kmax .or. kr.le.kmax)
               if (kl .gt. kmax) then
                  if (n(l-1) .ne. 0) then
                     lmax = l - 1
                  end if
                  do k=kr,kmax
                     m(l) = kl
                     n(l) = k
                     zl(l) = zi(i+1,j,k)
                     l = l + 1
                  end do
                  kr = kmax + 1
               else if (kr .gt. kmax) then
                  if (m(l-1) .ne. 0) then
                     lmax = l - 1
                  end if
                  do k=kl,kmax
                     m(l) = k
                     n(l) = kr
                     zl(l) = zi(i,j,k)
                     l = l + 1
                  end do
                  kl = kmax + 1
               else
                  m(l) = kl
                  n(l) = kr
                  if (zi(i,j,kl) .lt. zi(i+1,j,kr) ) then
                     zl(l) = zi(i,j,kl)
                     kl = kl + 1
!                    why not calc buoyr(kl) here?
                  else
                     zl(l) = zi(i+1,j,kr)
                     kr = kr + 1
!                    why not calc buoyl(kr) here?
                  end if
                  l = l + 1
               end if
            end do

            lmin = lextr + 1 ! should cause segmentation fault
            do l=1,lmax
               if (m(l).ne.0 .and. n(l).ne.0) then
                  lmin = l
                  exit
               end if
            end do            

            k = kmax
            l = lmax
            kl = m(lmax)
            kr = n(lmax)
            pdiff = _HALF_*(buoy(i,j,kmax)+buoy(i+1,j,kmax))*(ssen(i+1,j)-ssen(i,j))
            do while (k.ge.1 .and. l.ge.lmin)

               corr = _ZERO_
               zvel = _HALF_ * ( zz(i,j,k) + zz(i+1,j,k) )

               do while (l.ge.lmin .and. zl(l).gt.zvel)

                  zm = zz(i  ,j,m(l))
                  zn = zz(i+1,j,n(l))

!                 calculate dcm[l-1/2]
!                 calculate c[i+1](z[i,m(l)]) alias buoyplus
                  if (zm .ge. zz(i+1,j,kmax)) then
                     dcm = buoy(i+1,j,kmax)
                  else if (zm .le. zz(i+1,j,1)) then
                     dcm = buoy(i+1,j,1)
                  else
                     do while(zm .lt. zz(i+1,j,kr-1))
                        kr = kr - 1
                     end do
                     weight =   ( zm           - zz(i+1,j,kr-1) ) &
                              / ( zz(i+1,j,kr) - zz(i+1,j,kr-1) )
                     dcm =          weight * buoy(i+1,j,kr  ) &
                           + (_ONE_-weight)* buoy(i+1,j,kr-1)
                  end if
                  dcm = dcm - buoy(i,j,m(l))

!                 calculate dcn[l-1/2]
!                 calculate c[i](z[i+1,m(l)]) alias buoyminus
                  if (zn .ge. zz(i,j,kmax)) then
                     dcn = buoy(i,j,kmax)
                  else if (zn .le. zz(i,j,1)) then
                     dcn = buoy(i,j,1)
                  else
                     do while(zn .lt. zz(i,j,kl-1))
                        kl = kl - 1
                     end do
                     weight =   ( zn         - zz(i,j,kl-1) ) &
                              / ( zz(i,j,kl) - zz(i,j,kl-1) )
                     dcn =          weight * buoy(i,j,kl  ) &
                           + (_ONE_-weight)* buoy(i,j,kl-1)
                  end if
                  dcn = buoy(i+1,j,n(l)) - dcn

!                 calculate buoydiff
                  if (dcm*dcn .le. _ZERO_) then
                     buoydiff = _ZERO_
                  else if (dcm .gt. _ZERO_) then
                     buoydiff = min( dcm , dcn )
                  else
                     buoydiff = max( dcm , dcn )
                  end if

!                 calc p-gradient (even in the dzl of u-point)
                  pdiff = pdiff + (zl(l)-zl(l-1))*buoydiff

                  l = l - 1
                  if (zl(l) .lt. zvel) then
                     corr = ( zvel - zl(l) ) * buoydiff
                  end if
               end do

               idpdx(i,j,k) = _HALF_*(hn(i,j,k)+hn(i+1,j,k))*(pdiff-corr)*dxm1

               k = k - 1
            end do
            l = k
            do k=l,1,-1
               idpdx(i,j,k) = _ZERO_
            end do

         end if
      end do
   end do
!$OMP END DO


#ifndef SLICE_MODEL

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
            dym1=_ONE_/DYV
#endif
            l = 0
            kl = 0
            kr = 0
            lmax = -1 ! should cause segmentation fault
            do while(kl.le.kmax .or. kr.le.kmax)
               if (kl .gt. kmax) then
                  if (n(l-1) .ne. 0) then
                     lmax = l - 1
                  end if
                  do k=kr,kmax
                     m(l) = kl
                     n(l) = k
                     zl(l) = zi(i,j+1,k)
                     l = l + 1
                  end do
                  kr = kmax + 1
               else if (kr .gt. kmax) then
                  if (m(l-1) .ne. 0) then
                     lmax = l - 1
                  end if
                  do k=kl,kmax
                     m(l) = k
                     n(l) = kr
                     zl(l) = zi(i,j,k)
                     l = l + 1
                  end do
                  kl = kmax + 1
               else
                  m(l) = kl
                  n(l) = kr
                  if (zi(i,j,kl) .lt. zi(i,j+1,kr) ) then
                     zl(l) = zi(i,j,kl)
                     kl = kl + 1
                  else
                     zl(l) = zi(i,j+1,kr)
                     kr = kr + 1
                  end if
                  l = l + 1
               end if
            end do

            lmin = lextr+1 ! should cause segmentation fault
            do l=1,lmax
               if (m(l).ne.0 .and. n(l).ne.0) then
                  lmin = l
                  exit
               end if
            end do            

            k = kmax
            l = lmax
            kl = m(lmax)
            kr = n(lmax)
            pdiff = _HALF_*(buoy(i,j,kmax)+buoy(i,j+1,kmax))*(ssen(i,j+1)-ssen(i,j))
            do while (k.ge.1 .and. l.ge.lmin)

               corr = _ZERO_
               zvel = _HALF_ * ( zz(i,j,k) + zz(i,j+1,k) )

               do while (l.ge.lmin .and. zl(l).gt.zvel)

                  zm = zz(i,j  ,m(l))
                  zn = zz(i,j+1,n(l))

!                 calculate dcm[l-1/2]
!                 calculate c[i+1](z[i,m(l)]) alias buoyplus
                  if (zm .ge. zz(i,j+1,kmax)) then
                     dcm = buoy(i,j+1,kmax)
                  else if (zm .le. zz(i,j+1,1)) then
                     dcm = buoy(i,j+1,1)
                  else
                     do while(zm .lt. zz(i,j+1,kr-1))
                        kr = kr - 1
                     end do
                     weight =   ( zm           - zz(i,j+1,kr-1) ) &
                              / ( zz(i,j+1,kr) - zz(i,j+1,kr-1) )
                     dcm =          weight * buoy(i,j+1,kr  ) &
                           + (_ONE_-weight)* buoy(i,j+1,kr-1)
                  end if
                  dcm = dcm - buoy(i,j,m(l))

!                 calculate dcn[l-1/2]
!                 calculate c[i](z[i+1,m(l)]) alias buoyminus
                  if (zn .ge. zz(i,j,kmax)) then
                     dcn = buoy(i,j,kmax)
                  else if (zn .le. zz(i,j,1)) then
                     dcn = buoy(i,j,1)
                  else
                     do while(zn .lt. zz(i,j,kl-1))
                        kl = kl - 1
                     end do
                     weight =   ( zn         - zz(i,j,kl-1) ) &
                              / ( zz(i,j,kl) - zz(i,j,kl-1) )
                     dcn =          weight * buoy(i,j,kl  ) &
                           + (_ONE_-weight)* buoy(i,j,kl-1)
                  end if
                  dcn = buoy(i,j+1,n(l)) - dcn

!                 calculate buoydiff
                  if (dcm*dcn .le. _ZERO_) then
                     buoydiff = _ZERO_
                  else if (dcm .gt. _ZERO_) then
                     buoydiff = min( dcm , dcn )
                  else
                     buoydiff = max( dcm , dcn )
                  end if

!                 calc p-gradient (even in the dzl of v-point)
                  pdiff = pdiff + (zl(l)-zl(l-1))*buoydiff

                  l = l - 1
                  if (zl(l) .lt. zvel) then
                     corr = ( zvel - zl(l) ) * buoydiff
                  end if
               end do

               idpdy(i,j,k) = _HALF_*(hn(i,j,k)+hn(i,j+1,k))*(pdiff-corr)*dym1

               k = k - 1
            end do
            l = k
            do k=l,1,-1
               idpdy(i,j,k) = _ZERO_
            end do

         end if
      end do
   end do
!$OMP END DO

#endif


! Each thread must deallocate its own HEAP storage:
   deallocate(zl,stat=rc)
   if (rc /= 0) stop 'ip_stelling_vankester: Error deallocating memory (zl)'
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
