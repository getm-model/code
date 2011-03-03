!$Id: tow.F90,v 1.8 2009-05-25 19:00:18 hb Exp $
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
!  $Log: tow.F90,v $
!  Revision 1.8  2009-05-25 19:00:18  hb
!  correction for compiler option SLICE_MODEL
!
!  Revision 1.7  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.6  2007-05-26 12:19:30  kbk
!  destag of vertical velocities
!
!  Revision 1.5  2007-05-22 14:47:41  kbk
!  fixed indx for z-coordinates
!
!  Revision 1.4  2007-05-22 09:59:55  kbk
!  fixed HV index error + z-coordinate indices
!
!  Revision 1.3  2007-05-22 09:37:22  kbk
!  saving physical vertical velocities
!
!  Revision 1.2  2003-04-23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:20  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/04/21 09:41:33  bbh
!  Partial fixed problem with workspace (ws) variables in ncdf_save_?d.F90 and various conversion programs
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
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
