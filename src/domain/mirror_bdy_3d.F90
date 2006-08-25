!$Id: mirror_bdy_3d.F90,v 1.3 2006-08-25 09:00:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mirror_bdy_3d
!
! !INTERFACE:
   subroutine mirror_bdy_3d(f,tag)
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones, only : U_TAG,V_TAG,H_TAG,D_TAG
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
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
!  $Log: mirror_bdy_3d.F90,v $
!  Revision 1.3  2006-08-25 09:00:19  kbk
!  fixed sequence of boundary updates
!
!  Revision 1.2  2003/04/23 11:59:39  kbk
!  update_2d_halo on spherical variables + TABS to spaces
!
!  Revision 1.1  2003/04/07 15:22:03  kbk
!  parallel support
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
!
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
            f(nfi(n)-2,j,:) = f(nfi(n)-2,j-1,:)
            f(nfi(n)-1,j,:) = f(nfi(n)-1,j-1,:)
            do i = nfi(n),nli(n)
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j-1,:)
            end do
            f(nli(n)+1,j,:) = f(nfi(n)+1,j-1,:)
            f(nli(n)+2,j,:) = f(nfi(n)+2,j-1,:)
         end do
         do n = 1,NSB
            j = sj(n)
            f(sfi(n)-2,j,:) = f(sfi(n)-2,j+1,:)
            f(sfi(n)-1,j,:) = f(sfi(n)-1,j+1,:)
            do i = sfi(n),sli(n)
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j+1,:)
            end do
            f(sli(n)+1,j,:) = f(sli(n)+1,j+1,:)
            f(sli(n)+2,j,:) = f(sli(n)+2,j+1,:)
         end do
      case (V_TAG)
         do n = 1,NWB
            i = wi(n)
            f(i,wfj(n)-2,:) = f(i+1,wfj(n)-2,:)
            f(i,wfj(n)-1,:) = f(i+1,wfj(n)-1,:)
            do j = wfj(n),wlj(n)
               if (av(i,j) .eq. 3) f(i,j,:) = f(i+1,j,:)
            end do
            f(i,wlj(n)+1,:) = f(i+1,wlj(n)+1,:)
            f(i,wlj(n)+2,:) = f(i+1,wlj(n)+2,:)
         end do
         do n = 1,NEB
            i = ei(n)
            f(i,efj(n)-2,:) = f(i+1,efj(n)-2,:)
            f(i,efj(n)-1,:) = f(i+1,efj(n)-1,:)
            do j = efj(n),elj(n)
               if (av(i,j) .eq. 3) f(i,j,:) = f(i-1,j,:)
            end do
            f(i,elj(n)+1,:) = f(i+1,elj(n)+1,:)
            f(i,elj(n)+2,:) = f(i+1,elj(n)+2,:)
         end do
      case (H_TAG,D_TAG)
         do n = 1,NWB
            i = wi(n)
            f(i-1,wfj(n)-2,:) = f(i,wfj(n)-2,:)
            f(i-1,wfj(n)-1,:) = f(i,wfj(n)-1,:)
            do j = wfj(n),wlj(n)
               f(i-1,j,:) = f(i,j,:)
            end do
            f(i-1,wlj(n)+1,:) = f(i,wlj(n)+1,:)
            f(i-1,wlj(n)+2,:) = f(i,wlj(n)+2,:)
         end do
         do n = 1,NNB
            j = nj(n)
            f(nfi(n)-2,j+1,:) = f(nfi(n)-2,j,:)
            f(nfi(n)-1,j+1,:) = f(nfi(n)-1,j,:)
            do i = nfi(n),nli(n)
               f(i,j+1,:) = f(i,j,:)
            end do
            f(nli(n)+1,j+1,:) = f(nli(n)+1,j,:)
            f(nli(n)+2,j+1,:) = f(nli(n)+2,j,:)
         end do
         do n = 1,NEB
            i = ei(n)
            f(i+1,efj(n)-2,:) = f(i,efj(n)-2,:)
            f(i+1,efj(n)-1,:) = f(i,efj(n)-1,:)
            do j = efj(n),elj(n)
               f(i+1,j,:) = f(i,j,:)
            end do
            f(i+1,elj(n)+1,:) = f(i,elj(n)+1,:)
            f(i+1,elj(n)+2,:) = f(i,elj(n)+2,:)
         end do
         do n = 1,NSB
            j = sj(n)
            f(sfi(n)-2,j-1,:) = f(sfi(n)-2,j,:)
            f(sfi(n)-1,j-1,:) = f(sfi(n)-1,j,:)
            do i = sfi(n),sli(n)
               f(i,j-1,:) = f(i,j,:)
            end do
            f(sli(n)+1,j-1,:) = f(sli(n)+1,j,:)
            f(sli(n)+2,j-1,:) = f(sli(n)+2,j,:)
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
