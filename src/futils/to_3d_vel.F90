!$Id: to_3d_vel.F90,v 1.1.1.1 2002-05-02 14:01:20 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_3d_vel() - calculates 3D-velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_3d_vel(vel,trans,h,iimin,jjmin,kmin,iimax,jjmax,kmax,maxindx)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
  integer, intent(in)	:: maxindx
  integer, intent(in)	:: iimin,jjmin,kmin,iimax,jjmax,kmax
  REALTYPE, intent(in)	:: trans(I3DFIELD)
  REALTYPE, intent(in)	:: h(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)	:: vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: to_3d_vel.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:20  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/04/21 09:41:33  bbh
!  Partial fixed problem with workspace (ws) variables in ncdf_save_?d.F90 and various conversion programs
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   integer	:: indx
!EOP
!-----------------------------------------------------------------------
!BOC
   indx = 1
   do k=kmin,kmax
     do j=jjmin,jjmax
       do i=iimin,iimax
         vel(indx) = trans(i,j,k)/(h(i,j,k)+1.e-5)
         indx = indx+1
       end do
     end do
   end do

   return
   end subroutine to_3d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
