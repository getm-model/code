#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tow() - calculates vertical velocities and store in real*4.
!
! !INTERFACE:
   subroutine tow(imin,jmin,imax,jmax,kmin,kmax,mask,                  &
                  dt,                                                  &
#if defined CURVILINEAR || defined SPHERICAL
                  dxc,dyc,                                             &
#else
                  dx,dy,                                               &
#endif
                  HU,HV,hn,ho,uu,hun,vv,hvn,ww,missing,destag,ws)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: kmin(I2DFIELD)
   integer, intent(in)                 :: kmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: dt
#if defined CURVILINEAR || defined SPHERICAL
   REALTYPE, intent(in)                :: dxc(E2DFIELD)
   REALTYPE, intent(in)                :: dyc(E2DFIELD)
#else
   REALTYPE, intent(in)                :: dx,dy
#endif
   REALTYPE, intent(in)                :: HU(E2DFIELD)
   REALTYPE, intent(in)                :: HV(E2DFIELD)
   REALTYPE, intent(in)                :: hn(I3DFIELD)
   REALTYPE, intent(in)                :: ho(I3DFIELD)
   REALTYPE, intent(in)                :: uu(I3DFIELD)
   REALTYPE, intent(in)                :: hun(I3DFIELD)
   REALTYPE, intent(in)                :: vv(I3DFIELD)
   REALTYPE, intent(in)                :: hvn(I3DFIELD)
   REALTYPE, intent(in)                :: ww(I3DFIELD)
   REALTYPE, intent(in)                :: missing
   logical, intent(in)                 :: destag
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ws(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dtz,dxz,dyz
   REALTYPE                  :: u,v
   logical, save             :: physical_vel=.true.
!EOP
!-----------------------------------------------------------------------
!BOC

   if (physical_vel) then

!  save physical velocities
   do j=jmin,jmax
      do i=imin,imax
         if(mask(i,j) .gt. 0) then
!           bottom - normally k=0
            k=kmin(i,j)-1
            dtz=_ZERO_
            dxz=-(HU(i,j)-HU(i-1,j))/DXC
#ifndef SLICE_MODEL
            dyz=-(HV(i,j)-HV(i,j-1))/DYC
#else
            dyz=_ZERO_
#endif
            u=0.5*(uu(i,j,k+1)/hun(i,j,k+1)+uu(i-1,j,k+1)/hun(i-1,j,k+1))
            v=0.5*(vv(i,j,k+1)/hvn(i,j,k+1)+vv(i,j-1,k+1)/hvn(i,j-1,k+1))
            ws(i,j,k) = ww(i,j,k) + dtz + u*dxz + v*dyz
!           interior points
            do k=kmin(i,j),kmax-1
               dtz=dtz+(hn(i,j,k)-ho(i,j,k))/dt
               dxz=dxz+(hun(i,j,k)-hun(i-1,j,k))/DXC
#ifndef SLICE_MODEL
               dyz=dyz+(hvn(i,j,k)-hvn(i,j-1,k))/DYC
#else
               dyz=_ZERO_
#endif
               u=0.25*(uu(i,j,k  )/hun(i,j,k  )+uu(i-1,j,k  )/hun(i-1,j,k  )+&
                       uu(i,j,k+1)/hun(i,j,k+1)+uu(i-1,j,k+1)/hun(i-1,j,k+1) )
               v=0.25*(vv(i,j,k  )/hvn(i,j,k  )+vv(i,j-1,k  )/hvn(i,j-1,k  )+&
                       vv(i,j,k+1)/hvn(i,j,k+1)+vv(i,j-1,k+1)/hvn(i,j-1,k+1) )
               ws(i,j,k) = ww(i,j,k) + dtz + u*dxz + v*dyz
            end do
!           surface
            k=kmax
            dtz=dtz+(hn(i,j,k)-ho(i,j,k))/dt
            dxz=dxz+(hun(i,j,k)-hun(i-1,j,k))/DXC
#ifndef SLICE_MODEL
            dyz=dyz+(hvn(i,j,k)-hvn(i,j-1,k))/DYC
#else
               dyz=_ZERO_
#endif
            u=0.5*(uu(i,j,k)/hun(i,j,k)+uu(i-1,j,k)/hun(i-1,j,k))
            v=0.5*(vv(i,j,k)/hvn(i,j,k)+vv(i,j-1,k)/hvn(i,j-1,k))
            ws(i,j,k) = ww(i,j,k) + dtz  + u*dxz + v*dyz
            if (destag) then
               do k=kmax,kmin(i,j),-1
                  ws(i,j,k)=0.5*(ws(i,j,k)+ws(i,j,k-1))
               end do
               ws(i,j,k) = missing
            end if
         end if
      end do
   end do
   else
!     save grid-related velocities
      do k=0,kmax
         do j=jmin,jmax
            do i=imin,imax
               ws(i,j,k) = ww(i,j,k)
            end do
         end do
      end do
   end if

   return
   end subroutine tow
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
