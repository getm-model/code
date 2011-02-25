!$Id: eta_mask.F90,v 1.3 2003-05-09 11:38:26 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: eta_mask() - masks out land.
!
! !INTERFACE:
   subroutine eta_mask(imin,jmin,imax,jmax,mask,H,D,z,min_depth,missing, &
                       il,jl,ih,jh,eta)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: H(E2DFIELD)
   REALTYPE, intent(in)                :: D(E2DFIELD)
   REALTYPE, intent(in)                :: z(E2DFIELD)
   REALTYPE, intent(in)                :: min_depth,missing
   integer, intent(in)                 :: il,jl,ih,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: eta(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: eta_mask.F90,v $
!  Revision 1.3  2003-05-09 11:38:26  kbk
!  added proper undef support - based on Adolf Stips patch
!
!  Revision 1.2  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:18  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/04/21 09:41:33  bbh
!  Partial fixed problem with workspace (ws) variables in ncdf_save_?d.F90 and various conversion programs
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jl,jh
      do i=il,ih
         if (mask(i,j) .gt. 0 ) then
            if(z(i,j) .lt. (-H(i,j) + min_depth + 0.1) ) then
               eta(i,j) = missing
            else
               eta(i,j) = z(i,j)
            end if
         else
            eta(i,j) = missing
         end if
      end do
   end do

   return
   end subroutine eta_mask
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
