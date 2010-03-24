!$Id: divergence.F90,v 1.1 2010-03-24 14:58:15 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: divergence
!
! !INTERFACE:
   subroutine divergence()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
   use domain, only: ard1,dx,dy
#endif
   use variables_2d, only: surfdiv
   use variables_3d, only: hun,hvn,uu,vv
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: divergence.F90,v $
!  Revision 1.1  2010-03-24 14:58:15  kb
!  cleaned divergence calculations
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jmin,jmax
      do i=imin,imax
         surfdiv(i,j)= ((uu(i,  j  ,kmax)/hun(i  ,j  ,kmax)*DYU          &
                        -uu(i-1,j  ,kmax)/hun(i-1,j  ,kmax)*DYUIM1)      &
                       +(vv(i,  j  ,kmax)/hvn(i  ,j  ,kmax)*DXV          &
                        -vv(i,  j-1,kmax)/hvn(i  ,j-1,kmax)*DXVJM1))*ARCD1
      end do
   end do
   return
   end subroutine divergence
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
