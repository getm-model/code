#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion - 2D diffusion of momentum \label{sec-uv-diffusion}
!
! !INTERFACE:
   subroutine uv_diffusion(An_method,UEx,VEx,U,V,D,DU,DV, &
                           dudxC,dvdyC,shearX,AmC,AmX)

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
   use domain, only: dxc,dyc,dxx,dyx,arud1,arvd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_2d, only: PP,AnC,AnX
   use m2d, only: Am_method,Am_const,AM_CONSTANT,AM_LES,An_const
   use getm_timers,  only: tic,toc,TIM_UVDIFFUS
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                               :: An_method
   REALTYPE,dimension(E2DFIELD),intent(in),optional :: U,V,D,DU,DV
   REALTYPE,dimension(E2DFIELD),intent(in),optional :: dudxC,dvdyC,shearX
   REALTYPE,dimension(E2DFIELD),intent(in),optional :: AmC,AmX
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)       :: UEx,VEx
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!  Modified by       : Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE                                         :: depth
   integer                                          :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion() # ',Ncall
#endif
   CALL tic(TIM_UVDIFFUS)

#ifndef SLICE_MODEL
!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j)
#endif

!  Note (KK): see Kantha and Clayson 2000, page 591 for approximation
!             of diffusion terms
!             see Spurk 2008, page 484 for general orthogonal diffusion
!             see Griffies 2004, page 407 for transverse stress tensor
!             and resulting forces

!  Central for dx(2*Am*dx(U/DU))
#ifdef SLICE_MODEL
   j = jmax/2
#else
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
#endif
      do i=imin,imax+1          ! PP defined on T-points
         PP(i,j)=_ZERO_
!        Note (KK): we do not need N/S open boundary cells
!                   in W/E open boundary cells already dudxC=0
         if (az(i,j) .eq. 1) then
!           KK-TODO: center depth at velocity time stage
            select case(Am_method)
               case(AM_CONSTANT)
                  PP(i,j)=_TWO_*Am_const*DYC*D(i,j)*dudxC(i,j)
               case(AM_LES)
                  PP(i,j)=_TWO_*AmC(i,j)*DYC*D(i,j)*dudxC(i,j)
            end select
            select case(An_method)
!              KK-TODO: if gradient of velocity also accepted use of dudxC possible !?
               case(1)
                  PP(i,j)=PP(i,j)+An_const*DYC*(U(i,j)-U(i-1,j))/DXC
               case(2)
                  if (AnC(i,j) .gt. _ZERO_) then
                     PP(i,j)=PP(i,j)+AnC(i,j)*DYC*(U(i,j)-U(i-1,j))/DXC
                  end if
            end select
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#endif

#ifdef SLICE_MODEL
   j = jmax/2
#else
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
#endif
      do i=imin,imax      ! UEx defined on U-points
!        Note (KK): since U(au=3) will be obtained from mirroring
!                   we do not need PP in N/S open boundary cells
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j)=UEx(i,j)-(PP(i+1,j)-PP(i  ,j))*ARUD1
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   UEx(imin:imax,jmax/2+1) = UEx(imin:imax,jmax/2)
#endif

#ifndef SLICE_MODEL
!  Central for dy(Am*(dy(U/DU)+dx(V/DV)))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax        ! PP defined on X-points
         PP(i,j)=_ZERO_
         if (ax(i,j) .eq. 1) then
            select case(Am_method)
               case (AM_CONSTANT)
                  PP(i,j)=Am_const*DXX*_HALF_*(DU(i,j)+DU(i,j+1))*shearX(i,j)
               case (AM_LES)
                  PP(i,j)=AmX(i,j)*DXX*_HALF_*(DU(i,j)+DU(i,j+1))*shearX(i,j)
            end select
            select case(An_method)
!              Note (KK): outflow condition must be fulfilled !
!                         (use of mirrored transports or use of dudyX
!                          from deformation_rates, otherwise extended condition)
!                         at N/S closed boundaries slip condition dudyX=0
               case(1)
                  PP(i,j)=PP(i,j)+An_const*DXX*(U(i,j+1)-U(i,j))/DYX
               case(2)
                  if (AnX(i,j) .gt. _ZERO_) then
                     PP(i,j)=PP(i,j)+AnX(i,j)*DXX*(U(i,j+1)-U(i,j))/DYX
                  end if
            end select
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
        else if (av(i,j).eq.0 .and. av(i+1,j).eq.0) then
!           Note (KK): exclude convex corners (shearX=0)
!                      exclude W/E closed boundaries (not needed)
            if (au(i,j) .eq. 1) then ! northern closed boundary
               select case(Am_method)
                  case (AM_CONSTANT)
                     PP(i,j)=Am_const*DXX*DU(i,j)*shearX(i,j)
                  case (AM_LES)
                     PP(i,j)=AmX(i,j)*DXX*DU(i,j)*shearX(i,j)
               end select
            end if
            if (au(i,j+1) .eq. 1) then ! southern closed boundary
               select case(Am_method)
                  case (AM_CONSTANT)
                     PP(i,j)=Am_const*DXX*DU(i,j+1)*shearX(i,j)
                  case (AM_LES)
                     PP(i,j)=AmX(i,j)*DXX*DU(i,j+1)*shearX(i,j)
               end select
            end if
#endif
#endif
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax        !UEx defined on U-points
!        Note (KK): since U(au=3) will be obtained from mirroring
!                   we do not need PP at ax outside au=3
!                   at closed boundaries we need PP only at N/S closed boundaries
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j)=UEx(i,j)-(PP(i,j  )-PP(i,j-1))*ARUD1
         end if
      end do
   end do
!$OMP END DO
#endif

!  Central for dx(Am*(dy(U/DU)+dx(V/DV)))
#ifdef SLICE_MODEL
   j = jmax/2
#else
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
#endif
      do i=imin-1,imax      ! PP defined on X-points
         PP(i,j)=_ZERO_
         if (ax(i,j) .eq. 1) then
            select case(Am_method)
               case (AM_CONSTANT)
                  PP(i,j)=Am_const*DYX*_HALF_*(DV(i,j)+DV(i+1,j))*shearX(i,j)
               case (AM_LES)
                  PP(i,j)=AmX(i,j)*DYX*_HALF_*(DV(i,j)+DV(i+1,j))*shearX(i,j)
            end select
            select case(An_method)
!              Note (KK): outflow condition must be fulfilled !
!                         (use of mirrored transports or use of dvdxX
!                          from deformation_rates, otherwise extended condition)
!                         at W/E closed boundaries slip condition dvdxX=0
               case(1)
                  PP(i,j)=PP(i,j)+An_const*DYX*(V(i+1,j)-V(i,j))/DXX
               case(2)
                  if (AnX(i,j) .gt. _ZERO_) then
                     PP(i,j)=PP(i,j)+AnX(i,j)*DYX*(V(i+1,j)-V(i,j))/DXX
                  end if
            end select
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
        else if (au(i,j).eq.0 .and. au(i,j+1).eq.0) then
!           Note (KK): exclude convex corners (shearX=0)
!                      exclude N/S closed boundaries (not needed)
            if (av(i,j) .eq. 1) then ! eastern closed boundary
               select case(Am_method)
                  case (AM_CONSTANT)
                     PP(i,j)=Am_const*DXX*DV(i,j)*shearX(i,j)
                  case (AM_LES)
                     PP(i,j)=AmX(i,j)*DXX*DV(i,j)*shearX(i,j)
               end select
            end if
            if (av(i+1,j) .eq. 1) then ! western closed boundary
               select case(Am_method)
                  case (AM_CONSTANT)
                     PP(i,j)=Am_const*DXX*DV(i+1,j)*shearX(i,j)
                  case (AM_LES)
                     PP(i,j)=AmX(i,j)*DXX*DV(i+1,j)*shearX(i,j)
               end select
            end if
#endif
#endif
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#endif

#ifdef SLICE_MODEL
   j = jmax/2
#else
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
#endif
      do i=imin,imax          ! VEx defined on V-points
!        Note (KK): since V(av=3) will be obtained from mirroring
!                   we do not need PP at ax outside av=3
!                   at closed boundaries we need PP only at W/E closed boundaries
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j)=VEx(i,j)-(PP(i  ,j)-PP(i-1,j))*ARVD1
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   VEx(imin:imax,jmax/2-1) = VEx(imin:imax,jmax/2)
   VEx(imin:imax,jmax/2+1) = VEx(imin:imax,jmax/2)
#endif

#ifndef SLICE_MODEL
!  Central for dy(2*Am*dy(V/DV))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1     ! PP defined on T-points
      do i=imin,imax
         PP(i,j)=_ZERO_
!        Note (KK): we do not need W/E open boundary cells
!                   in N/S open boundary cells already dvdyC=0
         if (az(i,j) .eq. 1) then
!           KK-TODO: center depth at velocity time stage
            select case(Am_method)
               case (AM_CONSTANT)
                  PP(i,j)=_TWO_*Am_const*DXC*D(i,j)*dvdyC(i,j)
               case (AM_LES)
                  PP(i,j)=_TWO_*AmC(i,j)*DXC*D(i,j)*dvdyC(i,j)
            end select
            select case(An_method)
!              KK-TODO: if gradient of velocity also accepted use of dvdyC possible !?
               case(1)
                  PP(i,j)=PP(i,j)+An_const*DXC*(V(i,j)-V(i,j-1))/DYC
               case(2)
                  if (AnC(i,j) .gt. _ZERO_) then
                     PP(i,j)=PP(i,j)+AnC(i,j)*DXC*(V(i,j)-V(i,j-1))/DYC
                  end if
            end select
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax             ! VEx defined on V-points
      do i=imin,imax
!        Note (KK): since V(av=3) will be obtained from mirroring
!                   we do not need PP in W/E open boundary cells
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j)=VEx(i,j)-(PP(i,j+1)-PP(i,j  ))*ARVD1
         end if
      end do
   end do
!$OMP END DO
#endif

#ifndef SLICE_MODEL
!$OMP END PARALLEL
#endif

   CALL toc(TIM_UVDIFFUS)
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
