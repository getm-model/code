!$Id: divergence.F90,v 1.1.1.1 2002-05-02 14:00:43 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: divergence
!
! !INTERFACE:
   subroutine divergence()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use variables_2d, only: U,DU,V,DV,surfdiv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
      use domain, only: ard1,dx,dy
#endif
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
!  $Log: divergence.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:00:43  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.3  2001/05/25 18:48:12  bbh
!  Needs fixing - use mask
!
!  Revision 1.2  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j
!HB   REALTYPE	:: dxm1,dym1
!
!EOP
!-----------------------------------------------------------------------
!BOC
!HB   dxm1 = _ONE_/dx
!HB   dym1 = _ONE_/dy
return
   do i=imin,imax
      do j=jmin,jmax
      STDERR i,j,DU(i,j),DU(i-1,j),DV(i,j),DV(i,j-1)
         surfdiv(i,j) = ((U(i,j)/DU(i,j)*DYU-U(i-1,j)/DU(i-1,j)*DYUIM1)  &
                      - (V(i,j)/DV(i,j)*DXV-V(i,j-1)/DV(i,j-1)*DXVJM1))*ARCD1
      end do
   end do
   return
   end subroutine divergence
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
