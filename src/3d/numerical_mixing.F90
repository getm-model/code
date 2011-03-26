#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: numerical_mixing()
!
! !INTERFACE:
   subroutine numerical_mixing(F_2,F,nm3D,nm2d)
!
! !DESCRIPTION:
!
!
! !USES:
   use domain,       only: imin,imax,jmin,jmax,kmax
   use variables_3d, only: dt,hn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)  :: F_2(I3DFIELD)
   REALTYPE, intent(in)  :: F(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out) :: nm3d(I3DFIELD)
   REALTYPE, intent(out) :: nm2d(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   nm2d=_ZERO_
   do k=1,kmax
      do j=jmin,jmax
         do i=imin,imax
            nm3d(i,j,k)=(F_2(i,j,k)-F(i,j,k)**2)/dt
            nm2d(i,j)=nm2d(i,j)+nm3d(i,j,k)*hn(i,j,k)
         end do
      end do
   end do
   return
   end subroutine numerical_mixing
!BOC

!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hannes Rennau, Richard Hofmeister               !
!-----------------------------------------------------------------------
