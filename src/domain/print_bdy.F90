#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_bdy() - print open boundary info
!
! !INTERFACE:
   subroutine print_bdy(header)
!
! !DESCRIPTION:
!  Print the open boundary information. This routine is called twice -
!  first time with the global boundary infromation and second time
!  with the local boundary information. In the case of a serial run the
!  info is identical - in the case of a parallel run the open boundary 
!  information for a each sub-domain will be printed.
!
! !USES:
   use domain, only: NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: bdy_2d_type,bdy_3d_type
   use domain, only: bdy_2d_desc
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: header
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: m,n
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'print_bdy() # ',Ncall
#endif

   m=0
   LEVEL2 TRIM(header)
   LEVEL2 'There are ',NWB+NNB+NEB+NSB,' open boundaries.'
   if (NWB .ge. 1) then
      LEVEL3 'Western Boundary'
      do n = 1,NWB
         m=m+1
         LEVEL3 wi(n),wfj(n),wlj(n),trim(bdy_2d_desc(bdy_2d_type(m)))
      end do
   end if
   if (NNB .ge. 1) then
      LEVEL3 'Northern Boundary'
      do n = 1,NNB
         m=m+1
         LEVEL3 nj(n),nfi(n),nli(n),trim(bdy_2d_desc(bdy_2d_type(m)))
      end do
   end if
   if (NEB .ge. 1) then
      LEVEL3 'Eastern Boundary'
      do n = 1,NEB
         m=m+1
         LEVEL3 ei(n),efj(n),elj(n),trim(bdy_2d_desc(bdy_2d_type(m)))
      end do
   end if
   if (NSB .ge. 1) then
      LEVEL3 'Southern Boundary'
      do n = 1,NSB
         m=m+1
         LEVEL3 sj(n),sfi(n),sli(n),trim(bdy_2d_desc(bdy_2d_type(m)))
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving print_bdy()'
   write(debug,*)
#endif

   return
   end subroutine print_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
