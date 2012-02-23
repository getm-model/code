#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion - 2D diffusion of momentum \label{sec-uv-diffusion}
!
! !INTERFACE:
   subroutine uv_diffusion(An_method,UEx,VEx,U,V,D,DU,DV,hsd_u,hsd_v)

!  Note (KK): keep in sync with interface in m2d_general.F90
!
! !DESCRIPTION:
!
! Here, the diffusion terms for the vertically integrated transports are
! calculated by means of central differences, following the finite volume
! approach. They are added to the
! advection terms into the terms {\tt UEx} and {\tt VEx} for the
! $U$- and the $V$-equation, respectively. The physical diffusion with the
! given eddy viscosity coefficient $A_h^M$ is based on velocity gradients,
! whereas an additional numerical damping of the barotropic mode is based
! on gradients of the transports with the damping coefficient $A_h^N$,
! see the example given as equations (\ref{smooth_example_1}) and
! (\ref{smooth_example_2}).
!
! First diffusion term in (\ref{UMOM}):
! \begin{equation}
! \left(mn\,\partial_{\cal X}\left(2A_h^MD\partial_{\cal X}\left(\frac{U}{D}\right)
! + A_h^N\partial_{\cal X}U\right)\right)_{i,j}\approx
! \frac{
! {\cal F}^{Dxx}_{i+1,j}-{\cal F}^{Dxx}_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{equation}
!
! with diffusive fluxes
!
! \begin{equation}
! {\cal F}^{Dxx}_{i,j}=\left(2 A_h^MD_{i,j}\left(\frac{U_{i,j}}{D^u_{i,j}}
! -\frac{U_{i-1,j}}{D^u_{i-1,j}}\right)+A_h^N\left(U_{i,j}
! -U_{i-1,j}\right)\right)
! \frac{\Delta y^c_{i,j}}{\Delta x^c_{i,j}}.
! \end{equation}
!
! Second diffusion term in (\ref{UMOM}):
! \begin{equation}
! \left(mn\,\partial_{\cal Y}\left(A_h^MD\left(\partial_{\cal Y}\left(\frac{U}{D}\right)+\partial_{\cal X}\left(\frac{V}{D}\right)\right)
! + A_h^N\partial_{\cal Y}U\right)\right)_{i,j}\approx
! \frac{
! {\cal F}^{Dxy}_{i,j}-{\cal F}^{Dxy}_{i,j-1}
! }{\Delta x^x_{i,j}\Delta y^x_{i,j}}
! \end{equation}
!
! with diffusive fluxes
!
! \begin{equation}
! \begin{array}{rcl}
! \displaystyle
! {\cal F}^{Dxy}_{i,j}&=&
! \displaystyle
! A_h^M\frac12\left(D^u_{i,j}+D^u_{i,j+1}\right)
! \Delta x^x_{i,j}\left(\left(\frac{U_{i,j+1}}{D^u_{i,j+1}}
! -\frac{U_{i,j}}{D^u_{i,j}}\right)\frac{1}{\Delta y^x_{i,j}}
! +\left(\frac{V_{i+1,j}}{D^v_{i+1,j}}
! -\frac{V_{i,j}}{D^v_{i,j}}\right)\frac{1}{\Delta x^x_{i,j}}\right) \\ \\
! &&
! \displaystyle
! +A_h^N\left(U_{i,j+1} -U_{i,j}\right)\frac{\Delta x^x_{i,j}}
! {\Delta y^x_{i,j}}.
! \end{array}
! \end{equation}
!
! First diffusion term in (\ref{VMOM}):
! \begin{equation}
! \left(mn\,\partial_{\cal X}\left(A_h^MD\left(\partial_{\cal Y}\left(\frac{U}{D}\right)+\partial_{\cal X}\left(\frac{V}{D}\right)\right)
! + A_h^N\partial_{\cal X}V\right)\right)_{i,j}\approx
! \frac{
! {\cal F}^{Dyx}_{i,j}-{\cal F}^{Dyx}_{i-1,j}
! }{\Delta x^x_{i,j}\Delta y^x_{i,j}}
! \end{equation}
!
! with diffusive fluxes
!
! \begin{equation}
! \begin{array}{rcl}
! \displaystyle
! {\cal F}^{Dyx}_{i,j}&=&
! \displaystyle
! A_h^M\frac12\left(D^v_{i,j}+D^v_{i+1,j}\right)
! \Delta y^x_{i,j}\left(\left(\frac{U_{i,j+1}}{D^u_{i,j+1}}
! -\frac{U_{i,j}}{D^u_{i,j}}\right)\frac{1}{\Delta y^x_{i,j}}
! +\left(\frac{V_{i+1,j}}{D^v_{i+1,j}}
! -\frac{V_{i,j}}{D^v_{i,j}}\right)\frac{1}{\Delta x^x_{i,j}}\right) \\ \\
! &&
! \displaystyle
! +A_h^N\left(V_{i+1,j} -V_{i,j}\right)\frac{\Delta y^x_{i,j}}
! {\Delta x^x_{i,j}}.
! \end{array}
! \end{equation}
!
! Second diffusion term in (\ref{VMOM}):
! \begin{equation}
! \left(mn\,\partial_{\cal Y}\left(2A_h^MD\partial_{\cal Y}\left(\frac{V}{D}\right)
! + A_h^N\partial_{\cal Y}V\right)\right)_{i,j}\approx
! \frac{
! {\cal F}^{Dyy}_{i,j+1}-{\cal F}^{Dyy}_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{equation}
!
! with diffusive fluxes
!
! \begin{equation}
! {\cal F}^{Dyy}_{i,j}=\left(2 A_h^MD_{i,j}\left(\frac{V_{i,j}}{D^v_{i,j}}
! -\frac{V_{i,j-1}}{D^v_{i,j-1}}\right)+A_h^N\left(V_{i,j}
! -V_{i,j-1}\right)\right)
! \frac{\Delta x^c_{i,j}}{\Delta y^c_{i,j}}.
! \end{equation}
!
! The role of the additional diffusion of $U$ and $V$ with the
! diffusion coefficient $A_h^N$ is best demonstrated by means of a
! simplified set of vertically integrated equations:
!
! \begin{equation}\label{smooth_example_1}
! \begin{array}{l}
! \displaystyle
! \partial_t \eta = - \partial_x U - \partial_y V \\ \\
! \displaystyle
! \partial_t U = -gD\partial_x\eta
! + A_h^N \left(\partial_{xx} U + \partial_{yy} U\right) \\ \\
! \displaystyle
! \partial_t V = -gD\partial_y\eta
! + A_h^N \left(\partial_{xx} V + \partial_{yy} V\right), \\ \\
! \end{array}
! \end{equation}
!
! which can be transformed into an equation for $\partial_t\eta$ by
! derivation of the $\eta$-equation with respect to $t$,
! of the $U$-equation with respect to $x$
! and the $V$-equation with respect to $y$ and subsequent elimination of
! $U$ and $V$:
!
! \begin{equation}\label{smooth_example_2}
! \partial_t \left(\partial_t\eta\right) = gD \left(\partial_{xx}\eta
! + \partial_{yy}\eta\right)
! +  A_h^N \left(\partial_{xx}\left(\partial_t\eta\right)
! +\partial_{yy}\left(\partial_t\eta\right)\right),
! \end{equation}
!
! which can be interpreted as a wave equation with a damping on
! $\partial_t \eta$. This introduces an explicit damping of
! free surface elevation oscillations in a momentum-conservative
! manner. Hydrodynamic models with implicit treatment of the barotropic mode
! do not need to apply this method due to the implicit damping of those models,
! see e.g.\ \cite{BACKHAUS85}. The implementation of this explicit damping
! described here has been suggested by Jean-Marie Beckers, Li\'ege
! (Belgium).
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use m2d, only: Am
   use variables_2d, only: PP,An,AnX
   use getm_timers,  only: tic,toc,TIM_UVDIFF
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                                :: An_method
   REALTYPE,dimension(E2DFIELD),intent(in),optional  :: U,V,D,DU,DV
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)        :: UEx,VEx
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: hsd_u,hsd_v
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!  Modified by       : Knut Klingbeil
!
! !LOCAL VARIABLES:
   logical :: use_Am
   integer :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion() # ',Ncall
#endif
   CALL tic(TIM_UVDIFF)

   use_Am = (Am .gt. _ZERO_)

#ifdef _MOMENTUM_TERMS_
   if (present(hsd_u)) then
      hsd_u = UEx
   end if
   if (present(hsd_v)) then
      hsd_v = VEx
   end if
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j)

!  Central for dx(2*Am*dx(U/DU))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax+1          ! PP defined on T-points
         PP(i,j)=_ZERO_
!        Note (KK): we do not need N/S open bdy cells
!                   in W/E open bdy cells outflow condition must be
!                   explicitely prescribed (U(i,:) = U(i-1,:))
         if (az(i,j) .eq. 1) then
            if(use_Am) then
!              KK-TODO: we need center depth at velocity time stage
               PP(i,j)=_TWO_*Am*DYC*D(i,j)               &
                       *(U(i,j)/DU(i,j)-U(i-1,j)/DU(i-1,j))/DXC
            end if
            if (An_method .gt. 0) then
               PP(i,j)=PP(i,j)+An(i,j)*DYC*(U(i,j)-U(i-1,j))/DXC
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax      ! UEx defined on U-points
      do i=imin,imax
!        Note (KK): we do not need UEx(au=3)
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j)=UEx(i,j)-(PP(i+1,j)-PP(i  ,j))*ARUD1
         end if
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Central for dy(Am*(dy(U/DU)+dx(V/DV)))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax        ! PP defined on X-points
      do i=imin,imax
         PP(i,j)=_ZERO_
         if (ax(i,j) .ge. 1) then
            if(use_Am) then
               PP(i,j)=Am*_HALF_*(DU(i,j)+DU(i,j+1))*DXX  &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/DYX &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/DXX )
            end if
            if (An_method .gt. 0) then
               PP(i,j)=PP(i,j)+AnX(i,j)*(U(i,j+1)-U(i,j))*DXX/DYX
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax        !UEx defined on U-points
      do i=imin,imax
!        Note (KK): we do not need UEx(au=3)
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j)=UEx(i,j)-(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do
!$OMP END DO
#endif

!  Central for dx(Am*(dy(U/DU)+dx(V/DV)))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax      ! PP defined on X-points
      do i=imin-1,imax
         PP(i,j)=_ZERO_
         if (ax(i,j) .ge. 1) then
            if(use_Am) then
               PP(i,j)=Am*_HALF_*(DV(i,j)+DV(i+1,j))*DYX  &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/DYX &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/DXX )
            end if
            if (An_method .gt. 0) then
               PP(i,j)=PP(i,j)+AnX(i,j)*(V(i+1,j)-V(i,j))*DYX/DXX
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax          ! VEx defined on V-points
      do i=imin,imax
!        Note (KK): we do not need VEx(av=3)
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j)=VEx(i,j)-(PP(i  ,j)-PP(i-1,j))*ARVD1
         end if
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Central for dy(2*Am*dy(V/DV))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax
         PP(i,j)=_ZERO_
!        Note (KK): we do not need W/E open bdy cells
!                   in N/S open bdy cells outflow condition must be
!                   explicitely prescribed (V(:,j) = V(:,j-1))
         if (az(i,j) .eq. 1) then
            if(use_Am) then
!              KK-TODO: we need center depth at velocity time stage
               PP(i,j)=_TWO_*Am*DXC*D(i,j)               &
                       *(V(i,j)/DV(i,j)-V(i,j-1)/DV(i,j-1))/DYC
            end if
            if (An_method .gt. 0) then
               PP(i,j)=PP(i,j)+An(i,j)*DXC*(V(i,j)-V(i,j-1))/DYC
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax             ! VEx defined on V-points
      do i=imin,imax
!        Note (KK): we do not need VEx(av=3)
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j)=VEx(i,j)-(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do
!$OMP END DO
#endif

!$OMP END PARALLEL

#ifdef _MOMENTUM_TERMS_
   if (present(hsd_u)) then
      hsd_u = UEx - hsd_u
   end if
   if (present(hsd_v)) then
      hsd_v = VEx - hsd_v
   end if
#endif

   CALL toc(TIM_UVDIFF)
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
