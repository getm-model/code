#include "cppdefs.h"
!
   program test_speed_adv_upstream_2dh
!
! !DESCRIPTION:
!  program to test the eqstate inter action routines
!
!  To execute:
!  make test_speed_adv_upstream_2dh
!
! !USES:
   implicit none
!
! !LOCAL VARIABLES
!  basic eqstate interaction variables
   integer :: i,j,n
   integer, parameter :: Nmax=10000
   integer, parameter :: imin=1,imax=100,jmin=1,jmax=100

   REALTYPE :: Ah=_ZERO_
   REALTYPE :: dx=100.,dy=100.
   REALTYPE :: f(imin:imax,jmin:jmax)
   REALTYPE :: Dn(imin:imax,jmin:jmax)
   REALTYPE :: U(imin:imax,jmin:jmax)
   REALTYPE :: V(imin:imax,jmin:jmax)
   REALTYPE :: uflux(imin:imax,jmin:jmax)
   REALTYPE :: vflux(imin:imax,jmin:jmax)
   integer ::  az(imin:imax,jmin:jmax)
!EOP
!-----------------------------------------------------------------------
!BOC
   CALL RANDOM_NUMBER(f)

   do j=jmin,jmax
      do i=imin-1,imax
         if (f(i,j) .gt. 0.6666) then
            az(i,j) = 0
         else
            az(i,j) = 1
         end if
      end do
   end do

   do n=1,Nmax

   CALL RANDOM_NUMBER(f)
   CALL RANDOM_NUMBER(U)
   U = U - _HALF_
   CALL RANDOM_NUMBER(V)
   V = V - _HALF_

!#define _USE_ABS_
   
!  Calculating u-interface low-order fluxes !
   do j=jmin,jmax
      do i=imin-1,imax
#ifdef _USE_ABS_
         uflux(i,j) = _HALF_*( (U(i,j)+abs(U(i,j)))*f(i,j) + (U(i,j)-abs(U(i,j)))*f(i+1,j) ) 
#else
         if (U(i,j) .gt. _ZERO_) then
            uflux(i,j)=U(i,j)*f(i,j)
         else
            uflux(i,j)=U(i,j)*f(i+1,j)
         end if
#endif
#if 0
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0)) then
            uflux(i,j) = uflux(i,j) - AH*( f(i+1,j)                                 &
                                          -f(i  ,j))/DXU*_HALF_*(Dn(i+1,j)+Dn(i,j))
         end if
#endif
      end do
   end do

!  Calculating v-interface low-order fluxes !
   do j=jmin-1,jmax
      do i=imin,imax
#ifdef _USE_ABS_
         vflux(i,j) = _HALF_*( (V(i,j)+abs(V(i,j)))*f(i,j) + (V(i,j)-abs(V(i,j)))*f(i,j+1) ) 
#else
         if (V(i,j) .gt. _ZERO_) then
            vflux(i,j)=V(i,j)*f(i,j)
         else
            vflux(i,j)=V(i,j)*f(i,j+1)
         end if
#endif
#if 0
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0)) then
            vflux(i,j) = vflux(i,j) - AH*( f(i,j+1)                                 &
                                          -f(i,j  ))/DYV*_HALF_*(Dn(i,j+1)+Dn(i,j))
         end if
#endif
      end do
   end do

   end do

   end program test_speed_adv_upstream_2dh
!EOC

!-----------------------------------------------------------------------
! Copyright by the GETM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
