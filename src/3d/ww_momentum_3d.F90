!$Id: ww_momentum_3d.F90,v 1.2 2003-04-07 13:05:11 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ww_momentum_3d() - 3D-wmomentum equation.
!
! !INTERFACE:
   subroutine ww_momentum_3d()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,kmin,uu,vv,ww,ho,hn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: ww_momentum_3d.F90,v $
!  Revision 1.2  2003-04-07 13:05:11  kbk
!  cleaned code
!
!  Revision 1.1.1.1  2002/05/02 14:00:57  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.5  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.4  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.3  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.2  2001/05/01 07:13:27  bbh
!  use: kmax from m3d to domain
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   REALTYPE	:: dtm1
!HB   REALTYPE	:: dxm1,dym1
   integer	:: i,j,k
   REALTYPE	:: sold,snew,udiv,vdiv
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ww_momentum_3d() # ',Ncall
#endif

   dtm1=_ONE_/dt

   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (k .lt. kmin(i,j)) then
               ww(i,j,kmin(i,j))= _ZERO_
            else
               ww(i,j,k)=ww(i,j,k-1)			&
                        -(hn(i,j,k)-ho(i  ,j  ,k))*dtm1 &
                        -((uu(i,j,k)*DYU-uu(i-1,j  ,k)*DYUIM1) &
                        +(vv(i,j,k)*DXV-vv(i  ,j-1,k)*DXVJM1))*ARCD1
           end if
         end do
      end do
   end do

!     Consistency test: ww(i,j,kmax) must always be zero !

#ifdef DEBUG
   write(debug,*) 'Leaving ww_momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine ww_momentum_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
