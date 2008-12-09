!$Id: print_bdy.F90,v 1.3 2008-12-09 00:31:57 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_bdy(header) - prints boundary location info.
!
! !INTERFACE:
   subroutine print_bdy(header)
!
! !DESCRIPTION:
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: print_bdy.F90,v $
!  Revision 1.3  2008-12-09 00:31:57  kb
!  added new 2D open boundaries
!
!  Revision 1.2  2003-04-23 11:59:39  kbk
!  update_2d_halo on spherical variables + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:12  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: m,n
!
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
