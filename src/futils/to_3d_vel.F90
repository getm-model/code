!$Id: to_3d_vel.F90,v 1.5 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_vel() - calculates 3D-velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,mask, &
                        h,trans,missing,vel)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: kmin(I2DFIELD)
   integer, intent(in)                 :: kmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: h(I3DFIELD)
   REALTYPE, intent(in)                :: trans(I3DFIELD)
   REALTYPE, intent(in)                :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: vel(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: to_3d_vel.F90,v $
!  Revision 1.5  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.4  2006-01-11 14:02:42  lars
!  documentation + cosmetics
!
!  Revision 1.3  2003-05-09 11:38:26  kbk
!  added proper undef support - based on Adolf Stips patch
!
!  Revision 1.2  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:20  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/04/21 09:41:33  bbh
!  Partial fixed problem with workspace (ws) variables in ncdf_save_?d.F90 and various conversion programs
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=0,kmax
     do j=jmin,jmax
       do i=imin,imax
         if ( mask(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
            vel(i,j,k) = trans(i,j,k)/h(i,j,k)
         else
            vel(i,j,k) = missing
         end if
       end do
     end do
   end do
   return
   end subroutine to_3d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
