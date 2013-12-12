#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: filter_1d() - filter
!
! !INTERFACE:
   subroutine filter_1d(imin,imax,mask,var,il,ih,ws)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,imax
   integer, intent(in)                 :: mask(imin:imax)
   REALTYPE, intent(in)                :: var(imin:imax)
   integer, intent(in)                 :: il,ih
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: ws(imin:imax)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE,dimension(imin:imax) :: mean
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=imin,imax-1
      mean(i) = _HALF_*( var(i) + var(i+1) )
   end do

   do i=imin,il-1
      ws(i) = var(i)
   end do
   do i=il,ih
      if(mask(i) .eq. 1) then
         ws(i) = _HALF_*( mean(i-1) + mean(i) )
      else
         ws(i) = var(i)
      end if
   end do
   do i=ih+1,imax
      ws(i) = var(i)
   end do

   return
   end subroutine filter_1d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
