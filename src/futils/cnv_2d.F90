#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cnv_2d() - converts 2D-scalar fields to real*4.
!
! !INTERFACE:
   subroutine cnv_2d(imin,jmin,imax,jmax,mask,var,missing, &
                     il,jl,ih,jh,ws)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: var(E2DFIELD)
   REALTYPE, intent(in)                :: missing
   integer, intent(in)                 :: il,jl,ih,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: ws(E2DFIELD)
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
         if(mask(i,j) .gt. 0) then
            ws(i,j) = var(i,j)
         else
            ws(i,j) = missing
         end if
      end do
   end do
   return
   end subroutine cnv_2d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
