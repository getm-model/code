!$Id: depth_update.F90,v 1.3 2003-04-23 12:09:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: depth_update - adjust the depth to new elevations.
!
! !INTERFACE:
   subroutine depth_update
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,H,HU,HV,min_depth,crit_depth
   use domain, only: az,au,av,dry_z,dry_u,dry_v
   use variables_2d, only: D,z,zo,DU,zu,DV,zv
   use halo_zones, only : update_2d_halo,wait_halo,D_TAG,DU_TAG,DV_TAG
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
!  $Log: depth_update.F90,v $
!  Revision 1.3  2003-04-23 12:09:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 15:27:00  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:42  gotm
!  recovering after CVS crash
!
!  Revision 1.5  2001/08/27 11:53:13  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 13:03:34  bbh
!  Optimize for speed + masks for update_2d_halo() - CHECK
!
!  Revision 1.2  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: d1,d2,x
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'depth_update() # ',Ncall
#endif

#undef USE_MASK

!  U-points
   do j=jmin,jmax
      do i=imin,imax
#ifdef USE_MASK
         if(au(i,j) .gt. 0) then
#endif
         x=max(0.25*(zo(i,j)+zo(i+1,j)+z(i,j)+z(i+1,j)),-HU(i,j)+min_depth)
         zu(i,j) = x
         DU(i,j) = x+HU(i,j)
#ifdef USE_MASK
         end if
#endif
      end do
   end do

   call update_2d_halo(DU,DU,au,imin,jmin,imax,jmax,DU_TAG)

!  V-points
   do j=jmin,jmax
      do i=imin,imax
#ifdef USE_MASK
         if(av(i,j) .gt. 0) then
#endif
         x = max(0.25*(zo(i,j)+zo(i,j+1)+z(i,j)+z(i,j+1)),-HV(i,j)+min_depth)
         zv(i,j) = x
         DV(i,j) = x+HV(i,j)
#ifdef USE_MASK
         end if
#endif
      end do
   end do

   call wait_halo(DU_TAG)
!KBK   call cp_outside_openbdy_2d(DU)

   call update_2d_halo(DV,DV,av,imin,jmin,imax,jmax,DV_TAG)

!  Depth in elevation points
   D = z+H

   call wait_halo(DV_TAG)
!KBK   call cp_outside_openbdy_2d(DV)

   d1 = 2*min_depth
   d2 = crit_depth-2*min_depth
   where (az .gt. 0)
      dry_z = max(_ZERO_,min(_ONE_,(D-d1/2.)/d2))
   end where
   where (au .gt. 0)
      dry_u = max(_ZERO_,min(_ONE_,(DU-d1)/d2))
   end where
   where (av .gt. 0)
      dry_v = max(_ZERO_,min(_ONE_,(DV-d1)/d2))
   end where

#ifdef DEBUG
   do j=jmin,jmax
      do i=imin,imax

         if(D(i,j) .le. _ZERO_ .and. az(i,j) .gt. 0) then
            STDERR 'depth_update: D  ',i,j,H(i,j),D(i,j)
         end if

         if(DU(i,j) .le. _ZERO_ .and. au(i,j) .gt. 0) then
            STDERR 'depth_update: DU ',i,j,HU(i,j),DU(i,j)
         end if

         if(DV(i,j) .le. _ZERO_ .and. av(i,j) .gt. 0) then
            STDERR 'depth_update: DV ',i,j,HV(i,j),DV(i,j)
         end if

      end do
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving depth_update()'
   write(debug,*)
#endif
   return
   end subroutine depth_update
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
