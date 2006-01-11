!$Id: to_2d_vel.F90,v 1.4 2006-01-11 14:07:02 lars Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_2d_vel() - calculates 2D velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_2d_vel(imin,jmin,imax,jmax,mask,trans,depth,missing, &
                        il,jl,ih,jh,vel)
!
! !DESCRIPTION:
! This routine linearly interpolates the velocity at $u$-points to the $T$-points, 
! whenever the mask at the $T$-points is different from zero. Otherwise, the values
! are filled with the "missing value", {\tt missing}. The result is save in the 
! output argument {\tt vel}, which is single precision vector for storage in netCDF.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: trans(E2DFIELD)
   REALTYPE, intent(in)                :: depth(E2DFIELD)
   REALTYPE, intent(in)                :: missing
   integer, intent(in)                 :: il,jl,ih,jh
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)                ::  vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: to_2d_vel.F90,v $
!  Revision 1.4  2006-01-11 14:07:02  lars
!  documentation + cosmetics
!
!  Revision 1.3  2003-05-09 11:38:26  kbk
!  added proper undef support - based on Adolf Stips patch
!
!  Revision 1.2  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:19  gotm
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
   integer                   :: i,j,indx
!EOP
!-----------------------------------------------------------------------
!BOC
   indx = 1
   do j=jl,jh
      do i=il,ih
         if (mask(i,j) .gt. 0) then
            vel(indx) = trans(i,j)/depth(i,j)
         else
            vel(indx) = missing
         end if
         indx = indx+1
      end do
   end do

   return
   end subroutine to_2d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
