!$Id: slow_advection.F90,v 1.1 2002-05-02 14:00:54 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_advection() - ....
!
! !INTERFACE:
   subroutine slow_advection
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,HU,HV,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use m2d,    only: UEx,VEx,Uint,Vint
   use variables_3d,    only: ssun,ssvn
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
!  $Log: slow_advection.F90,v $
!  Revision 1.1  2002-05-02 14:00:54  gotm
!  Initial revision
!
!  Revision 1.6  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.5  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.4  2001/06/25 13:15:33  bbh
!  Fixed a few typos found by the DECFOR compiler
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,ii,jj
   REALTYPE	:: PP(iimin-1:iimax,jjmin-1:jjmax)
   REALTYPE	:: DUi(iimin-1:iimax+1,jjmin-1:jjmax+1)
   REALTYPE	:: DVi(iimin-1:iimax+1,jjmin-1:jjmax+1)
   REALTYPE	:: ps
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_advection() # ',Ncall
#endif

   do j=jjmin-1,jjmax+1
      do i=iimin-1,iimax+1
         DUi(i,j)=ssun(i,j)+HU(i,j)
      end do
   end do

   do j=jjmin-1,jjmax+1
      do i=iimin-1,iimax+1
         DVi(i,j)=ssvn(i,j)+HV(i,j)
      end do
   end do

! Upstream for dx(U^2/D)

   do j=jjmin,jjmax
      do i=iimin,iimax+1         ! PP defined on T-points
         if (az(i,j).ge.1) then
            PP(i,j)=0.5*(Uint(i-1,j)+Uint(i,j))
            if (PP(i,j) .gt. _ZERO_ ) then
	       ii=i-1
	    else
	       ii=i
	    end if
	    PP(i,j)=PP(i,j)*Uint(ii,j)/DUi(ii,j)*DYC
	 end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax           ! UEx defined on U-points
         if (au(i,j).ge.1) then
         UEx(i,j)=(PP(i+1,j)-PP(i  ,j))*ARUD1
	 end if
      end do
   end do

! Upstream for dy(UV/D)

   do i=iimin-1,iimax
      do j=jjmin-1,jjmax     ! PP defined on X-points
         if ((au(i,j).ge.1).or.(au(i,j+1).ge.1)) then
            PP(i,j)=0.5*(Vint(i+1,j)+Vint(i,j))
	    if (PP(i,j) .gt. _ZERO_) then
	       jj=j
	    else
	       jj=j+1
	    end if
	    PP(i,j)=PP(i,j)*Uint(i,jj)/DUi(i,jj)*DXX
	 end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax       !UEx defined on U-points
         if (au(i,j).ge.1) then
         UEx(i,j)=UEx(i,j)+(PP(i,j  )-PP(i,j-1))*ARUD1
	 end if
      end do
   end do

! Upstream for dx(UV/D)

   do j=jjmin-1,jjmax
      do i=iimin-1,iimax      ! PP defined on X-points
         if ((av(i,j).ge.1).or.(av(i+1,j).ge.1)) then
            PP(i,j)=0.5*(Uint(i,j)+Uint(i,j+1))
	    if (PP(i,j) .gt. _ZERO_) then
	       ii=i
	    else
	       ii=i+1
	    end if
	    PP(i,j)=PP(i,j)*Vint(ii,j)/DVi(ii,j)*DYX
	 end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax       ! VEx defined on V-points
         if (av(i,j).ge.1) then
         VEx(i,j)=(PP(i  ,j)-PP(i-1,j))*ARVD1
	 end if
      end do
   end do

! Upstream for dy(V^2/D)

   do j=jjmin,jjmax+1
      do i=iimin,iimax+1         ! PP defined on T-points
         if (az(i,j).ge.1) then
         PP(i,j)=0.5*(Vint(i,j-1)+Vint(i,j))
	 if (PP(i,j) .gt. _ZERO_) then
	    jj=j-1
	 else
	    jj=j
	 end if
	 PP(i,j)=PP(i,j)*Vint(i,jj)/DVi(i,jj)*DXC
	 end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax           ! VEx defined on V-points
         if (av(i,j).ge.1) then
         VEx(i,j)=VEx(i,j)+(PP(i,j+1)-PP(i,j  ))*ARVD1
	 end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving slow_advection()'
   write(debug,*)
#endif
   return
   end subroutine slow_advection
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
