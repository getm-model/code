!$Id: to_2d_vel.F90,v 1.2 2003-04-23 12:02:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_2d_vel() - calculates 2D velocities and store in real*4.
!
! !INTERFACE:
   subroutine to_2d_vel(vel,il,jl,ih,jh,mask,trans,depth,imin,jmin,imax,jmax)
!
! !DESCRIPTION:
!  Check version at JRC - KBK.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: il,jl,ih,jh
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: trans(E2DFIELD)
   REALTYPE, intent(in)                :: depth(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)::  vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: to_2d_vel.F90,v $
!  Revision 1.2  2003-04-23 12:02:43  kbk
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
         vel(indx) = trans(i,j)/(depth(i,j)+SMALL)
         indx = indx+1
      end do
   end do

   return
   end subroutine to_2d_vel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
