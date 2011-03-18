#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_vv() - average v-velocity to T-points
!
! !INTERFACE:
   subroutine to_3d_vv(imin,jmin,imax,jmax,kmin,kmax,az, &
                       hvn,vv,missing,vel)
!
! !DESCRIPTION:
! This routine linearly interpolates the velocity at $v$-points to the $T$-points, 
! whenever the mask at the $T$-points is different from zero. Otherwise, the values
! are filled with the "missing value", {\tt missing}. The result is written to the 
! output argument {\tt vel}, which is single precision vector for storage in netCDF.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)       :: imin,jmin,imax,jmax
  integer,  intent(in)       :: kmin(I2DFIELD)
  integer,  intent(in)       :: kmax
  integer,  intent(in)       :: az(E2DFIELD)
  REALTYPE, intent(in)       :: hvn(I3DFIELD)
  REALTYPE, intent(in)       :: vv(I3DFIELD)
  REALTYPE, intent(in)       :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)     :: vel(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: vt,vb
   REALTYPE, parameter       :: eps=1.E-5
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=0,kmax
      do j=jmin,jmax
         do i=imin,imax
         if ( az(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               vt        = vv(i,j  ,k)/(hvn(i,j  ,k)+eps)
               vb        = vv(i,j-1,k)/(hvn(i,j-1,k)+eps)
               vel(i,j,k) = 0.5*(vt+vb)
            else
               vel(i,j,k) = missing
            end if
         end do
      end do
   end do
   return
   end subroutine to_3d_vv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
