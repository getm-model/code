!$Id: uv_diffusion.F90,v 1.1 2002-05-02 14:00:47 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion() - diffusion of momentum.
!
! !INTERFACE:
   subroutine uv_diffusion(Am)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d,    only: D,U,DU,UEx,V,DV,VEx,PP
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  REALTYPE, intent(in) :: Am
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer 	:: i,j,ii,jj
   REALTYPE 	:: ps
   REALTYPE	:: p1,p2
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion() # ',Ncall
#endif

! Central for dx(2*Am*dx(uu^2/hun))
   do j=jmin,jmax
      do i=imin,imax+1          ! PP defined on T-points
         if (az(i,j) .ge. 1) then
            PP(i,j)=2.*Am*DYC*D(i,j)               &
               *(U(i,j)/DU(i,j)-U(i-1,j)/DU(i-1,j))/DXC
!HB            PP(i,j)=PP(i,j)+Am*DYC*(U(i,j)-U(i-1,j))/DXC
         end if
      end do
   end do

   do j=jmin,jmax      ! UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=UEx(i,j)-(PP(i+1,j)-PP(i  ,j))*ARUD1
	 end if
      end do
   end do

! Central for dy(Am*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jmin-1,jmax        ! PP defined on X-points
      do i=imin-1,imax
         if (au(i,j) .ge. 1 .or. au(i,j+1) .ge. 1) then
	    PP(i,j)=Am*0.5*(DU(i+1,j)+DU(i,j))*DXX  &
	          *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/DYX &
	           +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/DXX )
	 end if
      end do
   end do

   do j=jmin,jmax        !UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            UEx(i,j)=UEx(i,j)-(PP(i,j  )-PP(i,j-1))*ARUD1
	 end if
      end do
   end do

! Central for dx(Am*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jmin-1,jmax      ! PP defined on X-points
      do i=imin-1,imax
         if (av(i,j) .ge. 1 .or. av(i+1,j) .ge. 1) then
	    PP(i,j)=Am*0.5*(DU(i+1,j)+DU(i,j))*DXX  &
	       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/DYX &
	        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/DXX )
	 end if
      end do
   end do

   do j=jmin,jmax          ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=VEx(i,j)-(PP(i  ,j)-PP(i-1,j))*ARVD1
	 end if
      end do
   end do

! Central for dy(2*Am*dy(V^2/DV))
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax+1
	 if (az(i,j) .ge. 1) then
	    PP(i,j)=2.*Am*DXC*D(i,j)               &
	          *(V(i,j)/DV(i,j)-V(i,j-1)/DV(i,j-1))/DYC
!HB            PP(i,j)=PP(i,j)+Am*DXC*(V(i,j)-V(i,j-1))/DYC
	 end if
      end do
   end do

   do j=jmin,jmax             ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            VEx(i,j)=VEx(i,j)-(PP(i,j+1)-PP(i,j  ))*ARVD1
	 end if
      end do
   end do

#ifdef DEBUG
     write(debug,*) 'Leaving uv_diffusion()'
     write(debug,*)
#endif
   return
   end subroutine uv_diffusion
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
