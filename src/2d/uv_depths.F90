!$Id: uv_depths.F90,v 1.5 2003-04-23 12:09:44 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_depths - calculate depths in u and v points.
!
! !INTERFACE:
   subroutine uv_depths
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,H,HU,HV
   use variables_2d, only: DU,DV
   use halo_zones, only: update_2d_halo,wait_halo,HU_TAG,HV_TAG
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
!  $Log: uv_depths.F90,v $
!  Revision 1.5  2003-04-23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.4  2003/04/07 15:47:50  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:46  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.2  2001/05/18 12:55:13  bbh
!  Included masks in calls to update_2d_halo()
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
!  10Sep kbk: needs some more cleaning + wait for input from Hans
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   logical, save             :: first=.true.
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(0,*) 'uv_depths() # ',Ncall
#endif

   do j=jmin,jmax
      do i=imin,imax
#ifdef MIN_VEL_DEPTH
         HU(i,j)=min(H(i,j),H(i+1,j))
         HV(i,j)=min(H(i,j),H(i,j+1))
#else
         HU(i,j)=0.5*(H(i,j)+H(i+1,j))
         HV(i,j)=0.5*(H(i,j)+H(i,j+1))
#endif
      end do
   end do

   call update_2d_halo(HU,HU,au,imin,jmin,imax,jmax,HU_TAG)
   call wait_halo(HU_TAG)

   call update_2d_halo(HV,HV,av,imin,jmin,imax,jmax,HV_TAG)
   call wait_halo(HV_TAG)

#ifdef DEBUG
   write(0,*) 'Leaving uv_depths()'
   write(0,*)
#endif
   return
   end subroutine uv_depths
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
