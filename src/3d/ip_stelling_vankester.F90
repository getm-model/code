!$Id: ip_stelling_vankester.F90,v 1.1 2009-05-05 07:39:47 kb Exp $
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
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,l,kcount
   REALTYPE                  :: dxm1,dym1
   REALTYPE                  :: grdl,grdu,buoyl,prgr,dxz,dyz,dzz,zlm
   integer                   :: kplus,kminus,klower,kupper
   integer                   :: lnum=2*kmax+1
   REALTYPE                  :: zx(kmax)
   REALTYPE                  :: zl(2*kmax+1)
   REALTYPE                  :: dzl(2*kmax+1)
   REALTYPE                  :: dzfrac(kmax)
   REALTYPE                  :: db,dcn,dcm
   integer                   :: lvel(kmax+1)
   logical                   :: changed
   REALTYPE                  :: zltmp
   REALTYPE                  :: buoyplus,buoyminus
   REALTYPE                  :: zi(I3DFIELD)
   integer                   :: m(2*kmax+1)
   integer                   :: n(2*kmax+1)
!EOP
!-----------------------------------------------------------------------
!BOC
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!  First, the heights of the pressure points are calculated.
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+0.5*hn(i,j,1)
            zi(i,j,1)=-H(i,j)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+0.5*(hn(i,j,k-1)+hn(i,j,k))
               zi(i,j,k)=zi(i,j,k-1)+hn(i,j,k-1)
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
            ! merge zl vector from zi's
            zx(1)=-HU(i,j)+0.5*hun(i,j,1) ! zx defined on u-points
            zl(1)=zi(i,j,1)
            zl(2)=zi(i+1,j,1)
            do k=2,kmax
               zx(k)=zx(k-1)+0.5*(hun(i,j,k-1)+hun(i,j,k))
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

               if (zx(kcount) .gt. zl(l)) then
                  lvel(kcount)=l
                  dzfrac(kcount)=(zx(kcount)-zl(l))
                  kcount = kcount-1
               end if
            end do

            ! integrate internal pressure
            prgr=_ZERO_
            db=_ZERO_
            do k=kmax,1,-1
               do l=lvel(k+1)-1,lvel(k),-1
                  ! calc b-gradient
                  zlm=0.5*dzl(l)+zl(l)
                  
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

! Calculation of layer integrated internal pressure gradient as it
! appears on the right hand side of the v-velocity equation.
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dxm1 = _ONE_/DYV
#endif
            ! merge zl vector from zi's
            zx(1)=-HV(i,j)+0.5*hvn(i,j,1) ! zx defined on u-points
            zl(1)=zi(i,j,1)
            zl(2)=zi(i,j+1,1)
            do k=2,kmax
               zx(k)=zx(k-1)+0.5*(hvn(i,j,k-1)+hvn(i,j,k))
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

               if (zx(kcount) .gt. zl(l)) then
                  lvel(kcount)=l
                  dzfrac(kcount)=(zx(kcount)-zl(l))
                  kcount = kcount-1
               end if
            end do

            ! integrate internal pressure
            prgr=_ZERO_
            do k=kmax,1,-1
               do l=lvel(k+1)-1,lvel(k),-1
                  ! calc b-gradient
                  zlm=0.5*dzl(l)+zl(l)
                  
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

   return
   end subroutine ip_stelling_vankester
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
