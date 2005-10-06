!$Id: sealevel.F90,v 1.7 2005-10-06 09:54:00 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: sealevel() - using the cont. eq. to get the sealevel.
!
! !INTERFACE:
   subroutine sealevel
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,H
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only : arcd1,dxv,dyu
#else
   use domain, only : dx,dy,ard1
#endif
   use m2d, only: dtm
   use variables_2d, only: z,zo,U,V
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
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
!  $Log: sealevel.F90,v $
!  Revision 1.7  2005-10-06 09:54:00  hb
!  added support for vertical slice model - via -DSLICE_MODEL
!
!  Revision 1.6  2004/07/29 19:46:32  hb
!  For compiler option NOMADS_TEST: some lines shortened
!
!  Revision 1.5  2003/12/16 12:32:42  kbk
!  removed #ifdef SALTWEDGE_TEST (manuel)
!
!  Revision 1.4  2003/04/23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 15:44:13  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:45  gotm
!  recovering after CVS crash
!
!  Revision 1.9  2001/10/22 11:56:44  bbh
!  NOMADS_TEST kk =1.0 included
!
!  Revision 1.8  2001/09/01 17:15:13  bbh
!  Forgot to remove a few print statements
!
!  Revision 1.7  2001/09/01 16:49:27  bbh
!  A nasty hard-coding of z = 0. if mask = 2
!
!  Revision 1.6  2001/08/27 11:53:13  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.5  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 12:55:13  bbh
!  Included masks in calls to update_2d_halo()
!
!  Revision 1.2  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: kk
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'sealevel() # ',Ncall
#endif

   zo = z
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            z(i,j)=z(i,j)-dtm*((U(i,j)*DYU-U(i-1,j  )*DYUIM1) &
                              +(V(i,j)*DXV-V(i  ,j-1)*DXVJM1))*ARCD1

#ifdef NOMADS_TEST
       kk=1.0
       if ((((i.eq.1).or.(i.eq.imax)).and.(j.ge.1).and.(j.le.jmax)).or. &
           (((j.eq.1).or.(j.eq.jmax)).and.(i.ge.1).and.(i.le.imax)))    &
          z(i,j)=(1.-kk)*z(i,j)
       kk=0.5625
       if ((((i.eq.2).or.(i.eq.imax-1)).and.(j.ge.2).and.(j.le.jmax-1)).or. &
           (((j.eq.2).or.(j.eq.jmax-1)).and.(i.ge.2).and.(i.le.imax-1)))    &
          z(i,j)=(1.-kk)*z(i,j)
       kk=0.25
       if ((((i.eq.3).or.(i.eq.imax-2)).and.(j.ge.3).and.(j.le.jmax-2)).or. &
           (((j.eq.3).or.(j.eq.jmax-2)).and.(i.ge.3).and.(i.le.imax-2)))    &
           z(i,j)=(1.-kk)*z(i,j)
       kk=0.0625
       if ((((i.eq.4).or.(i.eq.imax-3)).and.(j.ge.4).and.(j.le.jmax-3)).or. &
           (((j.eq.4).or.(j.eq.jmax-3)).and.(i.ge.4).and.(i.le.imax-3)))    &
           z(i,j)=(1.-kk)*z(i,j)
#endif
         end if
      end do
   end do

#ifdef SLICE_MODEL
      do i=imin,imax
         z(i,3)=z(i,2)
      end do
#endif

   call update_2d_halo(z,z,az,imin,jmin,imax,jmax,z_TAG)
   call wait_halo(z_TAG)

#ifdef DEBUG
   write(debug,*) 'Leaving sealevel()'
   write(debug,*)
#endif
   return
   end subroutine sealevel
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
