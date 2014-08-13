#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: eta_mask() - masks out land.
!
! !INTERFACE:
   subroutine eta_mask(imin,jmin,imax,jmax,mask,H,D,z,mask_depth,missing, &
                       il,jl,ih,jh,eta)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: H(E2DFIELD)
   REALTYPE, intent(in)                :: D(E2DFIELD)
   REALTYPE, intent(in)                :: z(E2DFIELD)
   REALTYPE, intent(in)                :: mask_depth,missing
   integer, intent(in)                 :: il,jl,ih,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: eta(E2DFIELD)
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
         if (mask(i,j) .gt. 0 ) then
            if(z(i,j) .lt. (-H(i,j) + mask_depth) ) then
               eta(i,j) = missing
            else
               eta(i,j) = z(i,j)
            end if
         else
            eta(i,j) = missing
         end if
      end do
   end do

   return
   end subroutine eta_mask
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
