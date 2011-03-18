#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_vel() - calculates 3D-velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,mask, &
                        h,trans,missing,vel)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: kmin(I2DFIELD)
   integer, intent(in)                 :: kmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: h(I3DFIELD)
   REALTYPE, intent(in)                :: trans(I3DFIELD)
   REALTYPE, intent(in)                :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: vel(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=0,kmax
     do j=jmin,jmax
       do i=imin,imax
         if ( mask(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
            vel(i,j,k) = trans(i,j,k)/h(i,j,k)
         else
            vel(i,j,k) = missing
         end if
       end do
     end do
   end do
   return
   end subroutine to_3d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
