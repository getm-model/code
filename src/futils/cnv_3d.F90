!$Id: cnv_3d.F90,v 1.4 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cnv_3d() - convert 3D-scalar fields to real*4.
!
! !INTERFACE:
   subroutine cnv_3d(imin,jmin,imax,jmax,kmin,kmax,mask,var,missing, &
                     il,ih,jl,jh,kl,kh,ws)
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
   REALTYPE, intent(in)                :: var(I3DFIELD)
   REALTYPE, intent(in)                :: missing
   integer, intent(in)                 :: il,jl,ih,jh,kl,kh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ws(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: cnv_3d.F90,v $
!  Revision 1.4  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
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
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=kl,kh
      do j=jl,jh
         do i=il,ih
            if (mask(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               ws(i,j,k) = var(i,j,k)
            else
               ws(i,j,k) = missing
            end if
         end do
      end do
   end do
   return
   end subroutine cnv_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
