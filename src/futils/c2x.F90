#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Interpolate fields from central T-points to X-points
!
! !INTERFACE:
   subroutine c2x(imin,imax,jmin,jmax,cfield,xfield)
!
! !DESCRIPTION:
! This routine interpolates a variable given on the
! T-points to the X-points. At the edges and corners,
! an extrapolation is performed, since the outermost
! points are X-points in the GETM grid layout.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,    intent(in)              :: imin,imax,jmin,jmax
   REALTYPE,   intent(in)              :: cfield(I2DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE,  intent(out)              :: xfield(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: ufield(I2DFIELD),vfield(I2DFIELD)
 !EOP
!-----------------------------------------------------------------------
!BOC
!  do the interior X-points
   do j=jmin,jmax-1
      do i=imin,imax-1
         xfield(i,j) = 0.25*( cfield(i  ,j) + cfield(i+1,j+1) &
                           +  cfield(i+1,j) + cfield(i  ,j+1) )
      end do
   end do

!  do the interior U and V-points as an intermediate step
   do j=jmin,jmax
      do i=imin,imax-1
          ufield(i,j) = 0.5*( cfield(i,j) + cfield(i+1,j) )
      end do
   end do

   do j=jmin,jmax-1
      do i=imin,imax
         vfield(i,j) = 0.5*( cfield(i,j) + cfield(i,j+1) )
      end do
   end do

!  do the edges
   do i=imin,imax-1
      xfield(i,jmin-1) = 2.0*ufield(i,jmin) - xfield(i,jmin  )
      xfield(i,jmax  ) = 2.0*ufield(i,jmax) - xfield(i,jmax-1)
   end do

   do i=jmin,jmax-1
      xfield(imin-1,j) = 2.0*vfield(imin,j) - xfield(imin,j  )
      xfield(imax  ,j) = 2.0*vfield(imax,j) - xfield(imax-1,j)
   end do

!  do the exterior corners
   xfield(imin-1,jmin-1) = 2.0*ufield(imin-1,jmin) - xfield(imin-1,jmin  )
   xfield(imin-1,jmax  ) = 2.0*ufield(imin-1,jmax) - xfield(imin-1,jmax-1)

   xfield(imax  ,jmin-1) = 2.0*ufield(imax  ,jmin) - xfield(imax  ,jmin  )
   xfield(imax  ,jmax  ) = 2.0*ufield(imax  ,jmax) - xfield(imax  ,jmax-1)

   return
   end subroutine c2x

!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
