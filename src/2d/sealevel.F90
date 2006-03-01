!$Id: sealevel.F90,v 1.11 2006-03-01 15:54:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: sealevel - using the cont. eq. to get the sealevel.
!
! !INTERFACE:
   subroutine sealevel
!
! !DESCRIPTION:
!
! Here, the sea surface elevation is iterated according to the vertically
! integrated continuity equation given in (\ref{Elevation}) on page
! \pageref{Elevation}. 
!
! When working with the option {\tt SLICE\_MODEL}, the elevations
! at $j=2$ are copied to $j=3$.
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

#ifdef FRESHWATER_LENSE_TEST
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
