!$Id: cnv_3d.F90,v 1.1 2002-05-02 14:01:18 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cnv_3d() - convert 3D-scalar fields to real*4.
!
! !INTERFACE:
   subroutine cnv_3d(ws,var,iimin,jjmin,kmin,iimax,jjmax,kmax,maxindx)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: maxindx
   integer, intent(in)	:: iimin,jjmin,kmin,iimax,jjmax,kmax
   REALTYPE, intent(in)	:: var(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)	:: ws(*)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: cnv_3d.F90,v $
!  Revision 1.1  2002-05-02 14:01:18  gotm
!  Initial revision
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
            ws(indx) = var(i,j,k)
            indx = indx+1
         end do
      end do
   end do

   return
   end subroutine cnv_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
