!$Id: start_macro.F90,v 1.7 2004-01-05 12:43:22 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: start_macro() - initialize the macro loop.
!
! !INTERFACE:
   subroutine start_macro()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,H,HU,HV,min_depth
   use m2d, only: z,Uint,Vint
   use m3d, only: M
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
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
!  $Log: start_macro.F90,v $
!  Revision 1.7  2004-01-05 12:43:22  kbk
!  cleanned
!
!  Revision 1.6  2003/08/03 08:54:58  kbk
!  used HALO in loop boundaries
!
!  Revision 1.5  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.4  2003/04/07 16:27:32  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:55  gotm
!  recovering after CVS crash
!
!  Revision 1.5  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 13:00:36  bbh
!  Added masks for update_2d_halo() for possible later use
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: split
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'start_macro() # ',Ncall
#endif

   do j=jjmin-HALO,jjmax+HALO         ! Defining 'old' and 'new' sea surface
      do i=iimin-HALO,iimax+HALO      ! elevation for macro time step
         sseo(i,j)=ssen(i,j)
         ssen(i,j)=z(i,j)
      end do
   end do

   do j=jjmin-HALO,jjmax+HALO             ! Same for U-points
      do i=iimin-HALO,iimax+HALO-1
         ssuo(i,j)=ssun(i,j)
         ssun(i,j)=0.25*(sseo(i,j)+sseo(i+1,j)+ssen(i,j)+ssen(i+1,j))
         ssun(i,j)=max(ssun(i,j),-HU(i,j)+min_depth)
      end do
   end do

   do j=jjmin-HALO,jjmax+HALO-1
      do i=iimin-HALO,iimax+HALO             ! Same for V-points
         ssvo(i,j)=ssvn(i,j)
         ssvn(i,j)=0.25*(sseo(i,j)+sseo(i,j+1)+ssen(i,j)+ssen(i,j+1))
         ssvn(i,j)=max(ssvn(i,j),-HV(i,j)+min_depth)
      end do
   end do

! Defining vertically integrated, conservative
! u- and v-transport for macro time step

   split = _ONE_/float(M)
   Uint = split*Uint
   Vint = split*Vint

#ifdef DEBUG
   write(debug,*) 'Leaving start_macro()'
   write(debug,*)
#endif
   return
   end subroutine start_macro
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
