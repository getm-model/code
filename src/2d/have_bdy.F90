!$Id: have_bdy.F90,v 1.3 2003-04-23 12:09:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: have_bdy() - checks whether this node has boundaries.
!
! !INTERFACE:
   subroutine have_bdy
!
! !DESCRIPTION:
!
! !USES:
   use domain
   use m2d, only: have_boundaries
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: have_bdy.F90,v $
!  Revision 1.3  2003-04-23 12:09:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 15:42:05  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:44  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   integer                   :: nbdy
   integer                   :: f,l
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'have_bdy() # ',Ncall
#endif

   nbdy = 0
   i = 0
   if (NWB .ge. 1) then
      do n = 1,NWB
         if (wi(n) .ge. imin+ioff .and. wi(n) .le. imax+ioff) then
            wi(n) = wi(n) - ioff
            f = max(jmin+joff,wfj(n)) - joff
            l = min(jmax+joff,wlj(n)) - joff
            if(f .le. l) then
               i = i+1
               wfj(i) = f
               wlj(i) = l
               nbdy = nbdy+1
do k=1,nsbv
if(bdy_map(k,1) .eq. wi(n)+ioff .and. bdy_map(k,2) .eq. f+joff) then
bdy_index(nbdy) = k
end if
end do
            end if
         end if
      end do
   end if
   NWB = i

   i = 0
   if (NNB .ge. 1) then
      do n = 1,NNB
         if (nj(n) .ge. jmin+joff .and. nj(n) .le. jmax+joff) then
            nj(n) = nj(n) - joff
            f = max(imin+ioff,nfi(n)) - ioff
            l = min(imax+ioff,nli(n)) - ioff
            if(f .le. l) then
               i = i+1
               nfi(i) = f
               nli(i) = l
               nbdy = nbdy+1
do k=1,nsbv
if(bdy_map(k,1) .eq. f+ioff .and. bdy_map(k,2) .eq. nj(n)+joff) then
bdy_index(nbdy) = k
end if
end do
            end if
         end if
      end do
   end if
   NNB = i

   i = 0
   if (NEB .ge. 1) then
      do n = 1,NEB
         if (ei(n) .ge. imin+ioff .and. ei(n) .le. imax+ioff) then
            ei(n) = ei(n) - ioff
            f = max(jmin+joff,efj(n)) - joff
            l = min(jmax+joff,elj(n)) - joff
            if(f .le. l) then
               i = i+1
               efj(i) = f
               elj(i) = l
               nbdy = nbdy+1
do k=1,nsbv
if(bdy_map(k,1) .eq. ei(n)+ioff .and. bdy_map(k,2) .eq. l+joff) then
bdy_index(nbdy) = k
end if
end do
            end if
         end if
      end do
   end if
   NEB = i

   i = 0
   if (NSB .ge. 1) then
      do n = 1,NSB
         if (sj(n) .ge. jmin+joff .and. sj(n) .le. jmax+joff) then
            sj(n) = sj(n) - joff
            f = max(imin+ioff,sfi(n)) - ioff
            l = min(imax+ioff,sli(n)) - ioff
            if(f .le. l) then
               i = i+1
               sfi(i) = f
               sli(i) = l
               nbdy = nbdy+1
do k=1,nsbv
if(bdy_map(k,1) .eq. f+ioff .and. bdy_map(k,2) .eq. sj(n)+joff) then
bdy_index(nbdy) = k
end if
end do
            end if
         end if
      end do
   end if
   NSB = i

   if (nbdy .gt. 0) then
      have_boundaries = .true.
      bdy_index(nbdy+1:) = -1
   else
      have_boundaries = .false.
      bdy_index = -1
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving have_bdy()'
   write(debug,*)
#endif
   return
   end subroutine have_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
