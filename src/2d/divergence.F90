!$Id: divergence.F90,v 1.3 2003-04-23 12:09:43 kbk Exp $
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
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
      use domain, only: ard1,dx,dy
#endif
   use variables_2d, only: U,DU,V,DV,surfdiv
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
!  Revision 1.3  2003-04-23 12:09:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:08:37  kbk
!  cleaned code
!
!  Revision 1.1.1.1  2002/05/02 14:00:43  gotm
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
! !LOCAL VARIABLES:
   integer                   :: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=imin,imax
      do j=jmin,jmax
         surfdiv(i,j)=((U(i,j)/DU(i,j)*DYU-U(i-1,j)/DU(i-1,j)*DYUIM1)  &
                      -(V(i,j)/DV(i,j)*DXV-V(i,j-1)/DV(i,j-1)*DXVJM1)) &
                       *ARCD1
      end do
   end do
   return
   end subroutine divergence
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
