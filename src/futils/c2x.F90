!$Id: c2x.F90,v 1.1 2005-04-25 09:32:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Interpolate fields from central T-points to X-points
!
! !INTERFACE:
   subroutine c2x(iimin,iimax,jjmin,jjmax,cfield,xfield)
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
   integer,    intent(in)              :: iimin,jjmin,iimax,jjmax
   REALTYPE,   intent(in)              :: cfield(I2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE,  intent(out)              :: xfield(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: c2x.F90,v $
!  Revision 1.1  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: ufield(I2DFIELD),vfield(I2DFIELD)
 !EOP
!-----------------------------------------------------------------------
!BOC

!  do the interior X-points
   do j=jjmin,jjmax-1
      do i=iimin,iimax-1
         xfield(i,j) = 0.25*( cfield(i  ,j) + cfield(i+1,j+1) &
                           +  cfield(i+1,j) + cfield(i  ,j+1) )
      end do
   end do

!  do the interior U and V-points as an intermediate step
   do j=jjmin,jjmax
      do i=iimin,iimax-1
          ufield(i,j) = 0.5*( cfield(i,j) + cfield(i+1,j) ) 
      end do
   end do

   do j=jjmin,jjmax-1
      do i=iimin,iimax
         vfield(i,j) = 0.5*( cfield(i,j) + cfield(i,j+1) ) 
      end do
   end do

!  do the edges
   do i=iimin,iimax-1
      xfield(i,jjmin-1) = 2.0*ufield(i,jjmin) - xfield(i,jjmin  ) 
      xfield(i,jjmax  ) = 2.0*ufield(i,jjmax) - xfield(i,jjmax-1) 
   end do

   do i=jjmin,jjmax-1
      xfield(iimin-1,j) = 2.0*vfield(iimin,j) - xfield(iimin,j  ) 
      xfield(iimax  ,j) = 2.0*vfield(iimax,j) - xfield(iimax-1,j) 
   end do

!  do the exterior corners
   xfield(iimin-1,jjmin-1) = 2.0*ufield(iimin-1,jjmin) - xfield(iimin-1,jjmin  )
   xfield(iimin-1,jjmax  ) = 2.0*ufield(iimin-1,jjmax) - xfield(iimin-1,jjmax-1)

   xfield(iimax  ,jjmin-1) = 2.0*ufield(iimax  ,jjmin) - xfield(iimax  ,jjmin  )
   xfield(iimax  ,jjmax  ) = 2.0*ufield(iimax  ,jjmax) - xfield(iimax  ,jjmax-1)

   return
   end subroutine c2x

!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2005 - Lars Umlauf, Hans Burchard and Karsten Bolding
!-----------------------------------------------------------------------
