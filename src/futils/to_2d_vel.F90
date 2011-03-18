#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_2d_vel() - calculates 2D velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_2d_vel(imin,jmin,imax,jmax,mask,trans,depth,missing, &
                        il,jl,ih,jh,vel)
!
! !DESCRIPTION:
! This routine linearly interpolates the velocity at $u$-points to the $T$-points,
! whenever the mask at the $T$-points is different from zero. Otherwise, the values
! are filled with the "missing value", {\tt missing}. The result is save in the
! output argument {\tt vel}, which is single precision vector for storage in netCDF.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: trans(E2DFIELD)
   REALTYPE, intent(in)                :: depth(E2DFIELD)
   REALTYPE, intent(in)                :: missing
   integer, intent(in)                 :: il,jl,ih,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               ::  vel(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jl,jh
      do i=il,ih
         if (mask(i,j) .gt. 0) then
            vel(i,j) = trans(i,j)/depth(i,j)
         else
            vel(i,j) = missing
         end if
      end do
   end do
   return
   end subroutine to_2d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
