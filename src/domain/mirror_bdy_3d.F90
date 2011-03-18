#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mirror_bdy_3d() - mirrors 3d vaiables
!
! !INTERFACE:
   subroutine mirror_bdy_3d(f,tag)
!
! !DESCRIPTION:
!  Some variables are mirrored outside the calculation domain in the 
!  vicinity of the open boundaries. This is to avoid if statements
!  when calculating e.g. the Coriolis terms and advection.
!  This routines mirrors 3d variables.
!
! !USES:
   use halo_zones, only : U_TAG,V_TAG,H_TAG,D_TAG
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az,au,av
   use domain, only: NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: tag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'mirror_bdy_3d() # ',Ncall
#endif

   select case (tag)
      case (U_TAG)
         do n = 1,NNB
            j = nj(n)
            do i = nfi(n)-HALO,nli(n)+HALO
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j-1,:)
            end do
         end do
         do n = 1,NSB
            j = sj(n)
            do i = sfi(n)-HALO,sli(n)+HALO
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j+1,:)
            end do
         end do
      case (V_TAG)
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n)-HALO,wlj(n)+HALO
               if (av(i,j) .eq. 3) f(i,j,:) = f(i+1,j,:)
            end do
         end do
         do n = 1,NEB
            i = ei(n)
            do j = efj(n)-HALO,elj(n)+HALO
               if (av(i,j) .eq. 3) f(i,j,:) = f(i-1,j,:)
            end do
         end do
      case (H_TAG,D_TAG)
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n)-HALO,wlj(n)+HALO
               if (az(i,j) .gt. 1) f(i-1,j,:) = f(i,j,:)
            end do
         end do
         do n = 1,NNB
            j = nj(n)
            do i = nfi(n)-HALO,nli(n)+HALO
               if (az(i,j) .gt. 1) f(i,j+1,:) = f(i,j,:)
            end do
         end do
         do n = 1,NEB
            i = ei(n)
            do j = efj(n)-HALO,elj(n)+HALO
               if (az(i,j) .gt. 1) f(i+1,j,:) = f(i,j,:)
            end do
         end do
         do n = 1,NSB
            j = sj(n)
            do i = sfi(n)-HALO,sli(n)+HALO
               if (az(i,j) .gt. 1) f(i,j-1,:) = f(i,j,:)
            end do
         end do
      case default
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving mirror_bdy_3d()'
   write(debug,*)
#endif

   return
   end subroutine mirror_bdy_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2003 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
