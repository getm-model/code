#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_diffusion - slow diffusion terms \label{sec-slow-diffusion}
!
! !INTERFACE:
   subroutine slow_diffusion()
!
! !DESCRIPTION:
!
! Here, the calculation of the diffusive slow terms $S^x_D$ and $S^y_D$
! (see eqs.\ (\ref{SxD}) and (\ref{SxD})) is prepared. This routine
! basically repeats the calculations made in the routine {\tt uv\_iondvect},
! see section \ref{sec-uv-advect}, but this time based on the macro time
! step averaged and vertically integrated transports {\tt Uint} and {\tt Vint}.
! However, the damping of the external mode, as described in
! (\ref{smooth_example_1}) and
! (\ref{smooth_example_2}) is not considered here.
! The calculations of  $S^x_D$ and $S^y_D$ are then completed in the
! routine {\tt slow\_terms}, see section \ref{sec-slow-terms} on page
! \pageref{sec-slow-terms}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax,H,HU,HV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: D,U,V,UEx,VEx,Uint,Vint,PP
   use variables_3d, only: ssen,ssun,ssvn
#ifdef _LES_
   use variables_les, only: Am2d,AmX2d
#else
   use m2d, only: Am_const
#endif
   use getm_timers, only: tic, toc, TIM_SLOWDIFF
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii,jj
   REALTYPE                  :: Di(imin-1:imax+1,jmin-1:jmax+1)
   REALTYPE                  :: DUi(imin-1:imax+1,jmin-1:jmax+1)
   REALTYPE                  :: DVi(imin-1:imax+1,jmin-1:jmax+1)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_diffusion() # ',Ncall
#endif
   call tic(TIM_SLOWDIFF)

   do j=jmin-1,jmax+1
      do i=imin-1,imax+1
         Di(i,j)=ssen(i,j)+H(i,j)
      end do
   end do

   do j=jmin-1,jmax+1
      do i=imin-1,imax+1
         DUi(i,j)=ssun(i,j)+HU(i,j)
      end do
   end do

   do j=jmin-1,jmax+1
      do i=imin-1,imax+1
         DVi(i,j)=ssvn(i,j)+HV(i,j)
      end do
   end do

! Central for dx(2*AM*dx(U^2/HU))
   do j=jmin,jmax
      do i=imin,imax+1          ! PP defined on T-points
         PP(i,j)=_ZERO_
         if (az(i,j) .ge. 1) then
            if (_AM2D_ .gt. _ZERO_) then
               PP(i,j)=_TWO_*_AM2D_*DYC*Di(i,j)               &
                      *(Uint(i,j)/DUi(i,j)-Uint(i-1,j)/DUi(i-1,j))/DXC
            end if
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

#ifndef SLICE_MODEL
! Central for dy(AM*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jmin-1,jmax        ! PP defined on X-points
      do i=imin,imax
         PP(i,j)=_ZERO_
         if (ax(i,j) .ge. 1) then
            if (_AMX2D_ .gt. _ZERO_) then
               PP(i,j)=_AMX2D_*_HALF_*(DUi(i,j)+DUi(i,j+1))*DXX  &
                      *((Uint(i,j+1)/DUi(i,j+1)-Uint(i,j)/DUi(i,j))/DYX &
                       +(Vint(i+1,j)/DVi(i+1,j)-Vint(i,j)/DVi(i,j))/DXX )
            end if
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
#endif

! Central for dx(AM*(dy(U^2/DU)+dx(V^2/DV)))
   do j=jmin,jmax      ! PP defined on X-points
      do i=imin-1,imax
         PP(i,j)=_ZERO_
         if (ax(i,j) .ge. 1) then
            if (_AMX2D_ .gt. _ZERO_) then
               PP(i,j)=_AMX2D_*_HALF_*(DVi(i,j)+DVi(i+1,j))*DXX  &
                      *((Uint(i,j+1)/DUi(i,j+1)-Uint(i,j)/DUi(i,j))/DYX &
                       +(Vint(i+1,j)/DVi(i+1,j)-Vint(i,j)/DVi(i,j))/DXX )
            end if
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

#ifndef SLICE_MODEL
! Central for dy(2*AM*dy(V^2/DV))
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax
         PP(i,j)=_ZERO_
         if (az(i,j) .ge. 1) then
            if(_AM2D_ .gt. _ZERO_) then
               PP(i,j)=_TWO_*_AM2D_*DXC*Di(i,j)               &
                      *(Vint(i,j)/DVi(i,j)-Vint(i,j-1)/DVi(i,j-1))/DYC
            end if
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
#endif

   call toc(TIM_SLOWDIFF)
#ifdef DEBUG
     write(debug,*) 'Leaving slow_diffusion()'
     write(debug,*)
#endif
   return
   end subroutine slow_diffusion
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
