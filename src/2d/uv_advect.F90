!$Id: uv_advect.F90,v 1.4 2003-04-23 12:09:44 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_advect() - advection of momentum.
!
! !INTERFACE:
   subroutine uv_advect
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av
   use domain, only: ioff,joff
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: U,DU,UEx,V,DV,VEx,PP
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
!  $Log: uv_advect.F90,v $
!  Revision 1.4  2003-04-23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 15:58:18  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:46  gotm
!  recovering after CVS crash
!
!  Revision 1.7  2001/08/27 11:53:13  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.6  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.5  2001/06/25 13:15:33  bbh
!  Fixed a few typos found by the DECFOR compiler
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 13:03:34  bbh
!  Optimize for speed + masks for update_2d_halo() - CHECK
!
!  Revision 1.2  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii,jj
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect() # ',Ncall
#endif

!  Upstream for dx(U^2/D)
   do j=jmin,jmax         ! PP defined on T-points
      do i=imin,imax+1
         PP(i,j) = _ZERO_
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(U(i-1,j)+U(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               ii=i-1
            else
               ii=i
            end if
            PP(i,j)=PP(i,j)*U(ii,j)/DU(ii,j)*DYC
         end if
      end do
   end do
   do j=jmin,jmax      ! UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=(PP(i+1,j)-PP(i  ,j))*ARUD1
         else
            UEx(i,j)= _ZERO_
         end if
      end do
   end do

!  Upstream for dy(UV/D)
   do j=jmin-1,jmax        ! PP defined on X-points
      do i=imin-1,imax
         PP(i,j) = _ZERO_
         if (au(i,j) .ge. 1 .or. au(i,j+1) .ge. 1) then
            PP(i,j)=0.5*(V(i+1,j)+V(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j
            else
               jj=j+1
            end if
            PP(i,j)=PP(i,j)*U(i,jj)/DU(i,jj)*DXX
         end if
      end do
   end do
   do j=jmin,jmax        !UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=UEx(i,j)+(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do

! Upstream for dx(UV/D)
   do j=jmin-1,jmax      ! PP defined on X-points
      do i=imin-1,imax
         PP(i,j) = _ZERO_
         if (av(i,j) .ge. 1 .or. av(i+1,j) .ge. 1) then
            PP(i,j)=0.5*(U(i,j)+U(i,j+1))
            if (PP(i,j) .gt. _ZERO_) then
               ii=i
            else
               ii=i+1
            end if
            PP(i,j)=PP(i,j)*V(ii,j)/DV(ii,j)*DYX
         end if
      end do
   end do
   do j=jmin,jmax          ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=(PP(i  ,j)-PP(i-1,j))*ARVD1
         else
            VEx(i,j)= _ZERO_
         end if
      end do
   end do

!  Upstream for dy(V^2/D)
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax
         PP(i,j) = _ZERO_
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(V(i,j-1)+V(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j-1
            else
               jj=j
            end if
            PP(i,j)=PP(i,j)*V(i,jj)/DV(i,jj)*DXC
         end if
      end do
   end do
   do j=jmin,jmax             ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=VEx(i,j)+(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do

#ifdef DEBUG
     write(debug,*) 'Leaving uv_advect()'
     write(debug,*)
#endif
   return
   end subroutine uv_advect
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
