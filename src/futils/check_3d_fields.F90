#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: check_3d_fields() - check 3D-fields for having sane values
!
! !INTERFACE:
   subroutine check_3d_fields(imin,jmin,imax,jmax,kmin,kmax,mask,      &
                       f,minval,maxval,status)
!
! !DESCRIPTION:
!  This routine scans over a 3D field and checks that values are within
!  specified bounds. The routine is mainly intended for debugging
!  purposes - but can with little overhead also be applied for production
!  runs. If status is returned as a non-zero number out-of-bound values
!  have been found. It is up to the calling routine to take any actions.
!  The out-of-bound values are written to stderr.
!
! !USES:
   use getm_timers, only: tic, toc, TIM_CHECK3DF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)       :: imin,jmin,imax,jmax,kmax
  integer,  intent(in)       :: kmin(I2DFIELD)
  integer,  intent(in)       :: mask(E2DFIELD)
  REALTYPE, intent(in)       :: f(I3DFIELD)
  REALTYPE, intent(in)       :: minval,maxval
!
! !OUTPUT PARAMETERS:
   integer, intent(out)      :: status
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   call tic(TIM_CHECK3DF)
   status=0

   do k=0,kmax
      do j=jmin,jmax
         do i=imin,imax
            if ( mask(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               if(f(i,j,k) .lt. minval .or. f(i,j,k) .gt. maxval) then
                  LEVEL1 'check_3d_fields() ',i,j,k,f(i,j,k)
                  status=status+1
               end if
            end if
         end do
      end do
   end do
   call toc(TIM_CHECK3DF)

   return
   end subroutine check_3d_fields
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2006 - Karsten Bolding, Hans Burchard
!-----------------------------------------------------------------------
