!$Id: uv_advect.F90,v 1.12 2009-09-30 11:28:44 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_advect - 2D advection of momentum \label{sec-uv-advect}
!
! !INTERFACE:
   subroutine uv_advect
!
! !DESCRIPTION:
!
! The advective terms in the vertically integrated 
! momentum equation are discretised in
! a momentum-conservative form. This is carried out here for the 
! advective terms in the $U$-equation (\ref{UMOM}) and the 
! $V$-equation (\ref{VMOM}) (after applying the curvilinear
! coordinate transformationand multiplying these
! equations with $mn$). 
! 
! First advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{U^2}{Dn}\right)\right)_{i,j}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(U_{i+1,j}+U_{i,j})\tilde u_{i+1,j}\Delta y^c_{i+1,j}-
! \frac12(U_{i,j}+U_{i-1,j})\tilde u_{i,j}\Delta y^c_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
! 
! For the upwind scheme used here, the inter-facial velocities which are defined
! on T-points are here
! calculated as:
! 
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i-1,j}}{D^u_{i-1,j}} & \mbox{ for } \frac12(U_{i,j}+U_{i-1,j})>0\\ \\
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! Second advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle 
! \left(mn\,\partial_{\cal Y}y\left(\frac{UV}{Dm}\right)\right)_{i,j,k}\approx \\ \\ 
! \displaystyle 
! \quad
! \frac{
! \frac12(V_{i+1,j}+V_{i,j})\tilde u_{i,j}\Delta x^+_{i,j}-
! \frac12(V_{i+1,j-1}+V_{i,j-1})\tilde u_{i,j-1}\Delta x^+_{i,j-1}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
! 
! For the upwind scheme used here, the inter-facial 
! velocities which are defined on
! X-points are here
! calculated as:
! 
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ for } \frac12(V_{i+1,j,k}+V_{i,j,k})>0\\ \\
! \displaystyle
! \frac{U_{i,j+1}}{D^u_{i,j+1}} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! First advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle 
! \left(mn\,\partial_{\cal X}\left(\frac{UV}{Dn}\right)\right)_{i,j,k}\approx \\ \\ 
! \displaystyle 
! \quad
! \frac{
! \frac12(U_{i,j+1}+U_{i,j})\tilde v_{i,j}\Delta y^+_{i,j}-
! \frac12(U_{i-1,j+1}+U_{i-1,j})\tilde v_{i-1,j}\Delta y^+_{i-1,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
! 
! For the upwind scheme used here, the interfacial 
! velocities which are defined on
! X-points are here
! calculated as:
! 
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ for } \frac12(U_{i+1,j}+U_{i,j})>0\\ \\
! \displaystyle
! \frac{V_{i+1,j}}{D^v_{i+1,j}} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! Second advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}\left(\frac{V^2}{Dm}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(V_{i,j+1}+V_{i,j})\tilde v_{i,j+1}\Delta x^c_{i,j+1}-
! \frac12(V_{i,j}+V_{i,j-1})\tilde v_{i,j}\Delta x^c_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
! 
! For the upwind scheme used here, the interfacial velocities which are defined
! on T-points are here
! calculated as:
! 
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j-1}}{D^v_{i,j-1}} & \mbox{ for } \frac12(V_{i,j}+V_{i,j-1})>0\\ \\
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
   use domain, only: ioff,joff
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: U,DU,UEx,V,DV,VEx,PP
   use getm_timers, only: tic, toc, TIM_UVADVECT
!$ use omp_lib
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
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect() # ',Ncall
#endif
   CALL tic(TIM_UVADVECT)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj)

!  Upstream for dx(U^2/D)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax         ! PP defined on T-points
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(U(i-1,j)+U(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               ii=i-1
            else
               ii=i
            end if
            PP(i,j)=PP(i,j)*U(ii,j)/DU(ii,j)*DYC
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax      ! UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .eq. 1) then
            UEx(i,j)=(PP(i+1,j)-PP(i  ,j))*ARUD1
         end if
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Upstream for dy(UV/D)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax        ! PP defined on X-points
      do i=imin,imax
         if (ax(i,j) .ge. 1) then
            PP(i,j)=0.5*(V(i+1,j)+V(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j
            else
               jj=j+1
            end if
            PP(i,j)=PP(i,j)*U(i,jj)/DU(i,jj)*DXX
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax        !UEx defined on U-points
      do i=imin,imax
         if (au(i,j) .eq. 1) then
            UEx(i,j)=UEx(i,j)+(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do
!$OMP END DO
#endif

! Upstream for dx(UV/D)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax      ! PP defined on X-points
      do i=imin-1,imax
         if (ax(i,j) .ge. 1) then
            PP(i,j)=0.5*(U(i,j)+U(i,j+1))
            if (PP(i,j) .gt. _ZERO_) then
               ii=i
            else
               ii=i+1
            end if
            PP(i,j)=PP(i,j)*V(ii,j)/DV(ii,j)*DYX
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax          ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .eq. 1) then
            VEx(i,j)=(PP(i  ,j)-PP(i-1,j))*ARVD1
         end if
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Upstream for dy(V^2/D)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            PP(i,j)=0.5*(V(i,j-1)+V(i,j))
            if (PP(i,j) .gt. _ZERO_) then
               jj=j-1
            else
               jj=j
            end if
            PP(i,j)=PP(i,j)*V(i,jj)/DV(i,jj)*DXC
         else
            PP(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax             ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j) .eq. 1) then
            VEx(i,j)=VEx(i,j)+(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do
!$OMP END DO
#endif

!$OMP END PARALLEL

   CALL toc(TIM_UVADVECT)
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
