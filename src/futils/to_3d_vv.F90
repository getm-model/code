#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_vv() - average 3d-velocities to T-points
!
! !INTERFACE:
   subroutine to_3d_vv(vel,vv,hvn,av,iimin,jjmin,kmin,iimax,jjmax,kmax)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)       :: iimin,jjmin,kmin,iimax,jjmax,kmax
  REALTYPE, intent(in)       :: vv(I3DFIELD)
  REALTYPE, intent(in)       :: hvn(I3DFIELD)
  integer,  intent(in)       :: av(I2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out) :: vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: to_3d_vv.F90,v $
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   integer                   :: indx
   REALTYPE                  :: vb,vt
!EOP
!-----------------------------------------------------------------------
!BOC
   indx = 1
   do k=kmin,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            vb = _ZERO_
            vt = _ZERO_
            if (av(i,j-1) .ne. 0 .and. av(i,j) .ne. 0) then
               vb = vv(i,j-1,k)/(hvn(i,j-1,k)+1.e-5)
               vt = vv(i,j  ,k)/(hvn(i,j  ,k)+1.e-5)
            endif
            if (av(i,j-1) .eq. 0 .and. av(i,j) .ne. 0) then
               vt = vv(i,j  ,k)/(hvn(i,j  ,k)+1.e-5)
               vb = vt
            endif
            if (av(i,j-1) .ne. 0 .and. av(i,j) .eq. 0) then   
               vb = vv(i,j-1,k)/(hvn(i,j-1,k)+1.e-5)
               vt = vb
            endif
            vel(indx) = 0.5*(vb+vt)
            indx = indx+1
         end do
      end do
   end do

   return
   end subroutine to_3d_vv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
