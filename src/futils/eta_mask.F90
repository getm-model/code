!$Id: eta_mask.F90,v 1.1 2002-05-02 14:01:18 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: eta_mask() - masks out land.
!
! !INTERFACE:
   subroutine eta_mask(eta,il,jl,ih,jh,z,D,H,imin,jmin,imax,jmax,min_depth)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Check with version at JRC - KBK.
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: il,jl,ih,jh
   integer, intent(in)	:: imin,jmin,imax,jmax
   REALTYPE, intent(in)	:: z(E2DFIELD)
   REALTYPE, intent(in)	:: D(E2DFIELD)
   REALTYPE, intent(in)	:: H(E2DFIELD)
   REALTYPE, intent(in)	:: min_depth
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out) :: eta(*)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: eta_mask.F90,v $
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
   integer 	:: i,j
   integer 	:: indx
!EOP
!-----------------------------------------------------------------------
!BOC
   indx = 1
   do j=jl,jh
      do i=il,ih
         if(z(i,j) .lt. (-H(i,j) + min_depth + 0.1) ) then
            eta(indx) = -10.0
         else
            eta(indx) = z(i,j)
         end if
         indx = indx+1
      end do
   end do

   return
   end subroutine eta_mask
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
