#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_uu() - average u-velocities to T-points
!
! !INTERFACE:
   subroutine to_3d_uu(imin,jmin,imax,jmax,kmin,kmax,az, &
                       hun,uu,missing,vel)
!
! !DESCRIPTION:
! This routine linearly interpolates the velocity at $u$-points to the $T$-points, 
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
  REALTYPE, intent(in)       :: hun(I3DFIELD)
  REALTYPE, intent(in)       :: uu(I3DFIELD)
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
   REALTYPE                  :: ul,ur
   REALTYPE, parameter       :: eps=1.E-5
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=0,kmax
      do j=jmin,jmax
         do i=imin,imax
         if ( az(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               ul        = uu(i-1,j,k)/(hun(i-1,j,k)+eps)                     
               ur        = uu(i  ,j,k)/(hun(i  ,j,k)+eps)                     
               vel(i,j,k) = 0.5*(ul+ur)
            else
               vel(i,j,k) = missing
            end if
         end do
      end do
   end do
   return
   end subroutine to_3d_uu
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
