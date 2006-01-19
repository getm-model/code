#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_uu() - average u-velocities to T-points
!
! !INTERFACE:
   subroutine to_3d_uu(imin,jmin,imax,jmax,az,                          &
                       iimin,jjmin,iimax,jjmax,kmax,                    &
                       kmin,hun,uu,missing,vel)
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
  integer,  intent(in)       :: az(E2DFIELD)
  integer,  intent(in)       :: iimin,jjmin,iimax,jjmax,kmax
  integer,  intent(in)       :: kmin(I2DFIELD)
  REALTYPE, intent(in)       :: hun(I3DFIELD)
  REALTYPE, intent(in)       :: uu(I3DFIELD)
  REALTYPE, intent(in)       :: missing
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)      :: vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: to_3d_uu.F90,v $
!  Revision 1.4  2006-01-19 10:07:55  lars
!  bugfix (allocation of az-array)
!
!  Revision 1.3  2006-01-11 14:57:33  lars
!  added support for kmin
!
!  Revision 1.2  2006-01-11 14:03:28  lars
!  grave bug - partial re-write
!
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   integer                   :: indx
   REALTYPE                  :: ul,ur
   REALTYPE, parameter       :: eps=1.E-5
!EOP
!-----------------------------------------------------------------------
!BOC

   indx = 1
   do k=0,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
         if ( az(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               ul        = uu(i-1,j,k)/(hun(i-1,j,k)+eps)                     
               ur        = uu(i  ,j,k)/(hun(i  ,j,k)+eps)                     
               vel(indx) = 0.5*(ul+ur)
            else
               vel(indx) = missing
            end if
            indx = indx+1
         end do
      end do
   end do


   return
   end subroutine to_3d_uu
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
