!$Id: mirror_bdy_3d.F90,v 1.1 2003-04-07 15:22:03 kbk Exp $
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
   use halo_zones, only : U_TAG,V_TAG,H_TAG
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
!  Revision 1.1  2003-04-07 15:22:03  kbk
!  parallel support
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,n
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'mirror_bdy_3d() # ',Ncall
#endif
!KBK
   return
!KBK

   select case (tag)
      case (U_TAG)
#if 0
         do n = 1,NNB
            j = nj(n)
            do i = nfi(n),nli(n)
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j-1,:)
            end do
         end do
#else
STDERR 'mirror_bdy_3d: U_TAG'
#endif
         do n = 1,NSB
            j = sj(n)
            do i = sfi(n),sli(n)
               if (au(i,j) .eq. 3) f(i,j,:) = f(i,j+1,:)
            end do
         end do
      case (V_TAG)
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n),wlj(n)
               if (av(i,j) .eq. 3) f(i,j,:) = f(i+1,j,:)
            end do
         end do

         do n = 1,NEB
            i = ei(n)
            do j = efj(n),elj(n)
               if (av(i,j) .eq. 3) f(i,j,:) = f(i-1,j,:)
            end do
         end do
      case (H_TAG)
#if 1
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n),wlj(n)
               f(i-1,j,:) = f(i,j,:)
            end do
         end do

         do n = 1,NNB
            j = nj(n)
            do i = nfi(n),nli(n)
               f(i,j+1,:) = f(i,j,:)
            end do
         end do

         do n = 1,NEB
            i = ei(n)
            do j = efj(n),elj(n)
               f(i+1,j,:) = f(i,j,:)
            end do
         end do

         do n = 1,NSB
            j = sj(n)
            do i = sfi(n),sli(n)
               f(i,j-1,:) = f(i,j,:)
            end do
         end do
#endif
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
