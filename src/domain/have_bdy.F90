#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: have_bdy - checks whether this node has boundaries.
!
! !INTERFACE:
   subroutine have_bdy
!
! !DESCRIPTION:
!
! This routine which is called in {\tt domain.F90} checks whether the present
! node has open lateral boundaries. The integer field {\tt bdy\_index}
! is then set accordingly.
!
! !USES:
   use domain
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,m,n
   integer                   :: f,l
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'have_bdy() # ',Ncall
#endif

   i = 0
   m = 0
   if (NWB .ge. 1) then
      do n = 1,NWB
         m = m+1
         if (wi(n) .ge. imin+ioff .and. wi(n) .le. imax+ioff) then
            f = max(jlg,wfj(n)) - joff
            l = min(jhg,wlj(n)) - joff
            if(f.le.jmax .and. jmin.le.l) then
               nbdy = nbdy+1
               do k=bdy_index(m),nsbv
                  if (bdy_map(k,1) .eq. wi(n) .and. &
                      bdy_map(k,2) .eq. f+joff) then
                     bdy_index(nbdy) = k
                  end if
               end do
               bdy_2d_type(nbdy) = bdy_2d_type(m)
               bdy_3d_type(nbdy) = bdy_3d_type(m)
               bdy_index_l(nbdy) = nsbvl+1
               i = i+1
               wi(i) = wi(n) - ioff
               wfj(i) = f
               wlj(i) = l
               nsbvl = nsbvl + (l-f+1)
            end if
         end if
      end do
   end if
   NWB = i

   i = 0
   if (NNB .ge. 1) then
      do n = 1,NNB
         m = m+1
         if (nj(n) .ge. jmin+joff .and. nj(n) .le. jmax+joff) then
            f = max(ilg,nfi(n)) - ioff
            l = min(ihg,nli(n)) - ioff
            if(f.le.imax .and. imin.le.l) then
               nbdy = nbdy+1
               do k=bdy_index(m),nsbv
                  if (bdy_map(k,1) .eq. f+ioff .and.  &
                      bdy_map(k,2) .eq. nj(n)) then
                     bdy_index(nbdy) = k
                  end if
               end do
               bdy_2d_type(nbdy) = bdy_2d_type(m)
               bdy_3d_type(nbdy) = bdy_3d_type(m)
               bdy_index_l(nbdy) = nsbvl+1
               i = i+1
               nj(i) = nj(n) - joff
               nfi(i) = f
               nli(i) = l
               nsbvl = nsbvl + (l-f+1)
            end if
         end if
      end do
   end if
   NNB = i

   i = 0
   if (NEB .ge. 1) then
      do n = 1,NEB
         m = m+1
         if (ei(n) .ge. imin+ioff .and. ei(n) .le. imax+ioff) then
            f = max(jlg,efj(n)) - joff
            l = min(jhg,elj(n)) - joff
            if(f.le.jmax .and. jmin.le.l) then
               nbdy = nbdy+1
               do k=bdy_index(m),nsbv
                  if (bdy_map(k,1) .eq. ei(n) .and. &
                      bdy_map(k,2) .eq. f+joff) then
                     bdy_index(nbdy) = k
                  end if
               end do
               bdy_2d_type(nbdy) = bdy_2d_type(m)
               bdy_3d_type(nbdy) = bdy_3d_type(m)
               bdy_index_l(nbdy) = nsbvl+1
               i = i+1
               ei(i) = ei(n) - ioff
               efj(i) = f
               elj(i) = l
               nsbvl = nsbvl + (l-f+1)
            end if
         end if
      end do
   end if
   NEB = i

   i = 0
   if (NSB .ge. 1) then
      do n = 1,NSB
         m = m+1
         if (sj(n) .ge. jmin+joff .and. sj(n) .le. jmax+joff) then
            f = max(ilg,sfi(n)) - ioff
            l = min(ihg,sli(n)) - ioff
            if(f.le.imax .and. imin.le.l) then
               nbdy = nbdy+1
               do k=bdy_index(m),nsbv
                  if (bdy_map(k,1) .eq. f+ioff .and. &
                      bdy_map(k,2) .eq. sj(n)) then
                     bdy_index(nbdy) = k
                  end if
               end do
               bdy_2d_type(nbdy) = bdy_2d_type(m)
               bdy_3d_type(nbdy) = bdy_3d_type(m)
               bdy_index_l(nbdy) = nsbvl+1
               i = i+1
               sj(i) = sj(n) - joff
               sfi(i) = f
               sli(i) = l
               nsbvl = nsbvl + (l-f+1)
            end if
         end if
      end do
   end if
   NSB = i

   if (nbdy .gt. 0) then
      have_boundaries = .true.
      bdy_index(nbdy+1:) = -1
      bdy_index_l(nbdy+1:) = -1
   else
      have_boundaries = .false.
      bdy_index = -1
      bdy_index_l = -1
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
