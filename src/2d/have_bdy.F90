!$Id: have_bdy.F90,v 1.1.1.1 2002-05-02 14:00:44 gotm Exp $
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
!kbk   use m2d,    only: HaveBoundaries,EWBdy,ENBdy,EEBdy,ESBdy
   use m2d,    only: have_boundaries
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
!  Revision 1.1.1.1  2002-05-02 14:00:44  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer	:: i
   integer	:: n,nbdy
   integer	:: f,l
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
         if (wi(n) .ge. imin .and. wi(n) .le. imax) then
            if ((wfj(n) .ge. jmin .and. wfj(n) .le. jmax)  .or. &
                (wlj(n) .ge. jmin .and. wlj(n) .le. jmax)) then
               i = i+1
               f = max(jmin,wfj(n))
               wfj(i) = f
               l = min(jmax,wlj(n))
               wlj(i) = l
               nbdy = nbdy+1
            end if
         end if
      end do
   end if
   NWB = i

   i = 0
   if (NNB .ge. 1) then
      do n = 1,NNB
         if (nj(n) .ge. jmin .and. nj(n) .le. jmax) then
            if ((nfi(n) .ge. imin .and. nfi(n) .le. imax)  .or. &
                (nli(n) .ge. imin .and. nli(n) .le. imax)) then
               i = i+1
               f = max(imin,nfi(n))
               nfi(i) = f
               l = min(imax,nli(n))
               nli(i) = l
               nbdy = nbdy+1
            end if
         end if
      end do
   end if
   NNB = i

   i = 0
   if (NEB .ge. 1) then
      do n = 1,NEB
         if (ei(n) .ge. imin .and. ei(n) .le. imax) then
            if ((efj(n) .ge. jmin .and. efj(n) .le. jmax)  .or. &
                (elj(n) .ge. jmin .and. elj(n) .le. jmax)) then
               i = i+1
               f = max(jmin,efj(n))
               efj(i) = f
               l = min(jmax,elj(n))
               elj(i) = l
               nbdy = nbdy+1
            end if
         end if
      end do
   end if
   NEB = i

   i = 0
   if (NSB .ge. 1) then
      do n = 1,NSB
         if (sj(n) .ge. jmin .and. sj(n) .le. jmax) then
            if ((sfi(n) .ge. imin .and. sfi(n) .le. imax)   .or. &
                (sli(n) .ge. imin .and. sli(n) .le. imax)) then
               i = i+1
               f = max(imin,sfi(n))
               sfi(i) = f
               l = min(imax,sli(n))
               sli(i) = l
               nbdy = nbdy+1
            end if
         end if
      end do
   end if
   NSB = i

   if (nbdy .gt. 0) have_boundaries = .true.

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
