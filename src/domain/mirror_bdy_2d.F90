!$Id: mirror_bdy_2d.F90,v 1.5 2009-08-31 10:37:03 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mirror_bdy_2d
!
! !INTERFACE:
   subroutine mirror_bdy_2d(f,tag)
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones, only : U_TAG,V_TAG,H_TAG
   use domain, only: imin,imax,jmin,jmax
   use domain, only: az,au,av
   use domain, only: NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: tag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(E2DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: mirror_bdy_2d.F90,v $
!  Revision 1.5  2009-08-31 10:37:03  bjb
!  Consistent treatment of topo in halo zones
!
!  Revision 1.4  2007-05-14 08:12:43  kbk
!  fixed loops
!
!  Revision 1.3  2006-08-25 09:00:19  kbk
!  fixed sequence of boundary updates
!
!  Revision 1.2  2003/04/23 11:59:39  kbk
!  update_2d_halo on spherical variables + TABS to spaces
!
!  Revision 1.1  2003/04/07 15:22:03  kbk
!  parallel support
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'mirror_bdy_2d() # ',Ncall
#endif

   select case (tag)
      case (U_TAG)
         do n = 1,NNB
            j = nj(n)
            do i = nfi(n)-HALO,nli(n)+HALO
               if (au(i,j) .eq. 3) f(i,j) = f(i,j-1)
             end do
         end do
         do n = 1,NSB
            j = sj(n)
            do i = sfi(n)-HALO,sli(n)+HALO
               if (au(i,j) .eq. 3) f(i,j) = f(i,j+1)
            end do
         end do
      case (V_TAG)
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n)-HALO,wlj(n)+HALO
               if (av(i,j) .eq. 3) f(i,j) = f(i+1,j)
            end do
         end do
         do n = 1,NEB
            i = ei(n)
            do j = efj(n)-HALO,elj(n)+HALO
               if (av(i,j) .eq. 3) f(i,j) = f(i-1,j)
            end do
         end do
      case default
         do n = 1,NWB
            i = wi(n)
            do j = wfj(n)-HALO,wlj(n)+HALO
               if (az(i,j) .gt. 1) f(i-1,j) = f(i,j)
            end do
         end do

         do n = 1,NNB
            j = nj(n)
            do i = nfi(n)-HALO,nli(n)+HALO
               if (az(i,j) .gt. 1) f(i,j+1) = f(i,j)
            end do
         end do

         do n = 1,NEB
            i = ei(n)
            do j = efj(n)-HALO,elj(n)+HALO
               if (az(i,j) .gt. 1) f(i+1,j) = f(i,j)
            end do
         end do

         do n = 1,NSB
            j = sj(n)
            do i = sfi(n)-HALO,sli(n)+HALO
               if (az(i,j) .gt. 1) f(i,j-1) = f(i,j)
            end do
         end do
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving mirror_bdy_2d()'
   write(debug,*)
#endif

   return
   end subroutine mirror_bdy_2d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2003 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
