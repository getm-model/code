!$Id: slow_advection.F90,v 1.12 2009-08-18 10:24:44 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_advection - slow advection terms \label{sec-slow-advection}
!
! !INTERFACE:
   subroutine slow_advection
!
! !DESCRIPTION:
!
! Here, the calculation of the advective slow terms $S^x_A$ and $S^y_A$
! (see eqs.\ (\ref{SxA}) and (\ref{SyA})) is prepared. This routine
! basically repeats the calculations made in the routine {\tt uv\_advect}, 
! see section \ref{sec-uv-advect}, but this time based on the macro time
! step averaged and vertically integrated transports {\tt Uint} and {\tt Vint}.
! The calculations of  $S^x_A$ and $S^y_A$ are then completed in the 
! routine {\tt slow\_terms}, see section \ref{sec-slow-terms} on page 
! \pageref{sec-slow-terms}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,HU,HV,az,au,av,ax
   use domain, only: H,min_depth
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: UEx,VEx,Uint,Vint,PP
   use variables_3d, only: ssun,ssvn
   use getm_timers, only: tic, toc, TIM_SLOWADV
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
   integer                   :: i,j,ii,jj
   REALTYPE                  :: DUi(I2DFIELD)
   REALTYPE                  :: DVi(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_advection() # ',Ncall
#endif
   call tic(TIM_SLOWADV)

   do j=jmin-1,jmax+1
      do i=imin-1,imax+1
         DUi(i,j)=ssun(i,j)+HU(i,j)
         DVi(i,j)=ssvn(i,j)+HV(i,j)
      end do
   end do

! Upstream for dx(U^2/D)
   do j=jmin,jmax
      do i=imin,imax+1         ! PP defined on T-points
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(Uint(i-1,j)+Uint(i,j))
            if (PP(i,j) .gt. _ZERO_ ) then
               ii=i-1
            else
               ii=i
            end if
            PP(i,j)=PP(i,j)*Uint(ii,j)/DUi(ii,j)*DYC
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
   do j=jmin,jmax
      do i=imin,imax           ! UEx defined on U-points
         if (au(i,j) .eq. 1) then
            UEx(i,j)=(PP(i+1,j)-PP(i  ,j))*ARUD1
         end if
      end do
   end do

#ifndef SLICE_MODEL
!  Upstream for dy(UV/D)
   do j=jmin-1,jmax     ! PP defined on X-points
      do i=imin,imax
         if (ax(i,j) .ge. 1) then
            PP(i,j)=0.5*(Vint(i+1,j)+Vint(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j
            else
               jj=j+1
            end if
            PP(i,j)=PP(i,j)*Uint(i,jj)/DUi(i,jj)*DXX
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
   do j=jmin,jmax
      do i=imin,imax       !UEx defined on U-points
         if (au(i,j) .eq. 1) then
            UEx(i,j)=UEx(i,j)+(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do
#endif

! Upstream for dx(UV/D)
   do j=jmin,jmax
      do i=imin-1,imax      ! PP defined on X-points
         if (ax(i,j) .ge. 1) then
            PP(i,j)=0.5*(Uint(i,j)+Uint(i,j+1))
            if (PP(i,j) .gt. _ZERO_) then
               ii=i
            else
               ii=i+1
            end if
            PP(i,j)=PP(i,j)*Vint(ii,j)/DVi(ii,j)*DYX
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
   do j=jmin,jmax
      do i=imin,imax       ! VEx defined on V-points
         if (av(i,j) .eq. 1) then
            VEx(i,j)=(PP(i  ,j)-PP(i-1,j))*ARVD1
         end if
      end do
   end do

#ifndef SLICE_MODEL
!  Upstream for dy(V^2/D)
   do j=jmin,jmax+1          ! PP defined on T-points
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(Vint(i,j-1)+Vint(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j-1
            else
               jj=j
            end if
            PP(i,j)=PP(i,j)*Vint(i,jj)/DVi(i,jj)*DXC
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
   do j=jmin,jmax           ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .eq. 1) then
            VEx(i,j)=VEx(i,j)+(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do
#endif

   call toc(TIM_SLOWADV)
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
