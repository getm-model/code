!$Id: to_2d_u.F90,v 1.2 2006-01-11 14:02:42 lars Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Mask U-velocity and interpolate to T-points
!
! !INTERFACE:
   subroutine to_2d_u(imin,jmin,imax,jmax,az,u,DU,missing,              &
                      il,jl,ih,jh,vel)
!
! !DESCRIPTION:
! This routine linearly interpolates the vertically integrated velocity 
! at $U$-points to the $T$-points, whenever the mask at the $T$-points is different 
! from zero. Otherwise, the values are filled with a "missing value", {\tt missing}. 
! The result is written to the output argument {\tt vel}, which is single precision 
! vector for storage in netCDF.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: imin,jmin,imax,jmax
   integer,  intent(in)                :: az(E2DFIELD)
   REALTYPE, intent(in)                :: u(E2DFIELD)
   REALTYPE, intent(in)                :: DU(E2DFIELD)
   REALTYPE, intent(in)                :: missing
   integer,  intent(in)                :: il,jl,ih,jh

! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: vel(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: to_2d_u.F90,v $
!  Revision 1.2  2006-01-11 14:02:42  lars
!  documentation + cosmetics
!
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
! !LOCAL VARIABLES:
   integer                             :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jl,jh
      do i=il,ih
         if (az(i,j) .gt. 0) then
            vel(i,j) = 0.5*( u(i-1,j)/DU(i-1,j)                     &
                         +   u(i  ,j)/DU(i  ,j) )
         else
            vel(i,j) = missing
         end if
      end do
   end do
   return
   end subroutine to_2d_u
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
