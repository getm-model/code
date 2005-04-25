#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_uu() - average 3d-velocities to T-points
!
! !INTERFACE:
   subroutine to_3d_uu(vel,uu,hun,au,iimin,jjmin,kmin,iimax,jjmax,kmax)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)       :: iimin,jjmin,kmin,iimax,jjmax,kmax
  REALTYPE, intent(in)       :: uu(I3DFIELD)
  REALTYPE, intent(in)       :: hun(I3DFIELD)
  integer,  intent(in)       :: au(I2DFIELD)
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
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   integer                   :: indx
   REALTYPE                  :: ul,ur
!EOP
!-----------------------------------------------------------------------
!BOC

   indx = 1
   do k=kmin,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            ul = _ZERO_
            ur = _ZERO_
            if (au(i-1,j) .ne. 0 .and. au(i,j) .ne. 0) then
               ul = uu(i-1,j,k)/(hun(i-1,j,k)+1.e-5)
               ur = uu(i  ,j,k)/(hun(i  ,j,k)+1.e-5)
            endif
            if (au(i-1,j) .eq. 0 .and. au(i,j) .ne. 0) then
               ur = uu(i  ,j,k)/(hun(i  ,j,k)+1.e-5)
               ul = ur
            endif
            if (au(i-1,j) .ne. 0 .and. au(i,j) .eq. 0) then   
               ul = uu(i-1,j,k)/(hun(i-1,j,k)+1.e-5)
               ur = ul
            endif
            vel(indx) = 0.5*(ul+ur)
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
