!$Id: to_2d_v.F90,v 1.1 2005-04-25 09:32:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Mask V-velocity and interpolate to T-points
!
! !INTERFACE:
   subroutine to_2d_v(imin,jmin,imax,jmax,az,v,DV,missing,              &
                      il,jl,ih,jh,vel)
!
! !DESCRIPTION:
! This routine interpolates the $V$-velocity component from the $V$-points
! to the $T$-points. If the mask for the $T$-points is zero, the positions
! are filled up with "missing values". The result is written in a single-precision
! vector four netCDF output.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: imin,jmin,imax,jmax
   integer,  intent(in)                :: az(E2DFIELD)
   REALTYPE, intent(in)                :: v(E2DFIELD)
   REALTYPE, intent(in)                :: DV(E2DFIELD)
   REALTYPE, intent(in)                :: missing
   integer,  intent(in)                :: il,jl,ih,jh

! !OUTPUT PARAMETERS:
   REAL_4B, intent(out)                :: vel(*)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: to_2d_v.F90,v $
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                             :: i,j,indx
!EOP
!-----------------------------------------------------------------------
!BOC
   indx = 1
   do j=jl,jh
      do i=il,ih
         if (az(i,j) .gt. 0) then
            vel(indx) = 0.5*( v(i,j-1)/DV(i,j-1)                        &
                          +   v(i,j  )/DV(i,j  ) )
         else
            vel(indx) = missing
         end if
         indx = indx+1
      end do
   end do

   return
   end subroutine to_2d_v
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------