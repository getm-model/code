#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_diff_2dh - lateral diffusion of velocity
! \label{sec-uv-diff-2dh}
!
! !INTERFACE:
   subroutine uv_diff_2dh(An_method,UEx,VEx,U,V,D,DU,DV,  &
                          dudxC,dvdyC,dudyX,dvdxX,shearX, &
                          AmC,AmX,phydis,hsd_u,hsd_v)
!  Note (KK): keep in sync with interface in m2d.F90
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
   use variables_2d, only: AnC,AnX
   use m2d, only: Am_method,Am_const,NO_AM,AM_LAPLACE,AM_LES,AM_CONSTANT,An_const
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                                :: An_method
   REALTYPE,dimension(E2DFIELD),intent(in),optional  :: U,V,D,DU,DV
   REALTYPE,dimension(E2DFIELD),intent(in),optional  :: dudxC,dvdyC
   REALTYPE,dimension(:,:),pointer,intent(in),optional  :: dudyX,dvdxX
   REALTYPE,dimension(E2DFIELD),intent(in),optional  :: shearX
   REALTYPE,dimension(E2DFIELD),intent(in),optional  :: AmC,AmX
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)        :: UEx,VEx
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(:,:),pointer,intent(out)       :: phydis
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: hsd_u,hsd_v
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!  Modified by       : Knut Klingbeil
!
! !LOCAL VARIABLES:
   logical                      :: calc_phydis
   REALTYPE,dimension(E2DFIELD) :: work2d
   REALTYPE,dimension(E2DFIELD) :: phydis_wrk,phydis_vel
   integer                      :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diff_2dh() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif

   calc_phydis = associated(phydis)

   if (present(hsd_u)) then
      hsd_u = UEx
   end if
   if (present(hsd_v)) then
      hsd_v = VEx
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)

!  Note (KK): see Kantha and Clayson 2000, page 591 for approximation
!             of diffusion terms
!             see Spurk 2008, page 484 for general orthogonal diffusion
!             see Griffies 2004, page 407 for transverse stress tensor
!             and resulting forces


!  diffusion terms for u-equation

   if (calc_phydis) then
      phydis_vel = _ZERO_
   end if

!  Central for dx(2*Am*dx(U/DU))
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin,jmax
#endif
      do i=imin-1,imax+1          ! work2d defined on T-points
         work2d(i,j)=_ZERO_
         if (az(i,j) .eq. 1) then
            select case(Am_method)
               case(AM_LAPLACE)
                  work2d(i,j)=Am_const*dudxC(i,j)
               case(AM_LES)
                  work2d(i,j)=_TWO_*AmC(i,j)*dudxC(i,j)
               case(AM_CONSTANT)
                  work2d(i,j)=_TWO_*Am_const*dudxC(i,j)
            end select
            if (Am_method .ne. NO_AM) then
               if (calc_phydis) then
                  phydis_wrk(i,j) = work2d(i,j) * dudxC(i,j)
               end if
               work2d(i,j) = DYC * D(i,j) * work2d(i,j)
            end if
            select case(An_method)
               case(1)
                  work2d(i,j)=work2d(i,j)+An_const*DYC*(U(i,j)-U(i-1,j))/DXC
               case(2)
                  if (AnC(i,j) .gt. _ZERO_) then
                     work2d(i,j)=work2d(i,j)+AnC(i,j)*DYC*(U(i,j)-U(i-1,j))/DXC
                  end if
            end select
         else
            if (calc_phydis) then
               phydis_wrk(i,j) = _ZERO_
            end if
         end if
      end do
#ifdef SLICE_MODEL
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
#endif
      do i=imin-1,imax      ! UEx defined on U-points
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j)=UEx(i,j)-(work2d(i+1,j)-work2d(i  ,j))*ARUD1
            if (calc_phydis) then
               phydis_vel(i,j) = _HALF_*(phydis_wrk(i,j)+phydis_wrk(i+1,j))
            end if
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

#ifndef SLICE_MODEL
   if (An_method.gt.0 .or. Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
!     Central for dy(Am*(dy(U/DU)+dx(V/DV)))
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax        ! work2d defined on X-points
         do i=imin-1,imax
            work2d(i,j)=_ZERO_
            if (calc_phydis) then
               phydis_wrk(i,j) = _ZERO_
            end if
            if (ax(i,j) .ge. 1) then
               select case(Am_method)
                  case (AM_CONSTANT)
                     work2d(i,j)=Am_const*shearX(i,j)
                  case (AM_LES)
                     work2d(i,j)=AmX(i,j)*shearX(i,j)
               end select
               if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                  if (calc_phydis) then
                     phydis_wrk(i,j) = work2d(i,j) * dudyX(i,j)
                  end if
                  work2d(i,j) = DXX * _HALF_*(DU(i,j)+DU(i,j+1)) * work2d(i,j)
               end if
               select case(An_method)
!                 Note (KK): outflow condition must be fulfilled !
!                            (use of mirrored transports or use of dudyX
!                             from deformation_rates, otherwise extended condition)
!                            at N/S closed boundaries slip condition dudyX=0
                  case(1)
                     work2d(i,j)=work2d(i,j)+An_const*DXX*(U(i,j+1)-U(i,j))/DYX
                  case(2)
                     if (AnX(i,j) .gt. _ZERO_) then
                        work2d(i,j)=work2d(i,j)+AnX(i,j)*DXX*(U(i,j+1)-U(i,j))/DYX
                     end if
               end select
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
            else if (av(i,j).eq.0 .and. av(i+1,j).eq.0) then
!              Note (KK): exclude convex corners (shearX=0)
!                         exclude W/E closed boundaries (not needed)
               if (au(i,j) .eq. 1) then ! northern closed boundary
                  select case(Am_method)
                     case (AM_CONSTANT)
                        work2d(i,j)=Am_const*shearX(i,j)
                     case (AM_LES)
                        work2d(i,j)=AmX(i,j)*shearX(i,j)
                  end select
                  if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                     if (calc_phydis) then
                        phydis_wrk(i,j) = work2d(i,j) * dudyX(i,j)
                     end if
                     work2d(i,j) = DXX * DU(i,j) * work2d(i,j)
                  end if
               end if
               if (au(i,j+1) .eq. 1) then ! southern closed boundary
                  select case(Am_method)
                     case (AM_CONSTANT)
                        work2d(i,j)=Am_const*shearX(i,j)
                     case (AM_LES)
                        work2d(i,j)=AmX(i,j)*shearX(i,j)
                  end select
                  if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                     if (calc_phydis) then
                        phydis_wrk(i,j) = work2d(i,j) * dudyX(i,j)
                     end if
                     work2d(i,j) = DXX * DU(i,j+1) * work2d(i,j)
                  end if
               end if
#endif
#endif
            end if
         end do
      end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax        !UEx defined on U-points
         do i=imin-1,imax
            if (au(i,j).eq.1 .or. au(i,j).eq.2) then
               UEx(i,j)=UEx(i,j)-(work2d(i,j  )-work2d(i,j-1))*ARUD1
               if (calc_phydis) then
                  phydis_vel(i,j) = phydis_vel(i,j) + _HALF_*(phydis_wrk(i,j-1)+phydis_wrk(i,j))
               end if
            end if
         end do
      end do
!$OMP END DO
   end if
#endif

   if (calc_phydis) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               phydis(i,j) = _HALF_*(phydis_vel(i-1,j)+phydis_vel(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
   end if


!  diffusion terms for v-equation

   if (calc_phydis) then
      phydis_vel = _ZERO_
   end if

   if (An_method.gt.0 .or. Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
!     Central for dx(Am*(dy(U/DU)+dx(V/DV)))
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-1,jmax
#endif
         do i=imin-1,imax      ! work2d defined on X-points
            work2d(i,j)=_ZERO_
            if (calc_phydis) then
               phydis_wrk(i,j) = _ZERO_
            end if
            if (ax(i,j) .ge. 1) then
               select case(Am_method)
                  case (AM_CONSTANT)
                     work2d(i,j)=Am_const*shearX(i,j)
                  case (AM_LES)
                     work2d(i,j)=AmX(i,j)*shearX(i,j)
               end select
               if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                  if (calc_phydis) then
                     phydis_wrk(i,j) = work2d(i,j) * dvdxX(i,j)
                  end if
                  work2d(i,j) = DYX * _HALF_*(DV(i,j)+DV(i+1,j)) * work2d(i,j)
               end if
               select case(An_method)
!                 Note (KK): outflow condition must be fulfilled !
!                            (use of mirrored transports or use of dvdxX
!                             from deformation_rates, otherwise extended condition)
!                            at W/E closed boundaries slip condition dvdxX=0
                  case(1)
                     work2d(i,j)=work2d(i,j)+An_const*DYX*(V(i+1,j)-V(i,j))/DXX
                  case(2)
                     if (AnX(i,j) .gt. _ZERO_) then
                        work2d(i,j)=work2d(i,j)+AnX(i,j)*DYX*(V(i+1,j)-V(i,j))/DXX
                     end if
               end select
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
           else if (au(i,j).eq.0 .and. au(i,j+1).eq.0) then
!              Note (KK): exclude convex corners (shearX=0)
!                         exclude N/S closed boundaries (not needed)
               if (av(i,j) .eq. 1) then ! eastern closed boundary
                  select case(Am_method)
                     case (AM_CONSTANT)
                        work2d(i,j)=Am_const*shearX(i,j)
                     case (AM_LES)
                        work2d(i,j)=AmX(i,j)*shearX(i,j)
                  end select
                  if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                     if (calc_phydis) then
                        phydis_wrk(i,j) = work2d(i,j) * dvdxX(i,j)
                     end if
                     work2d(i,j) = DXX * DV(i,j) * work2d(i,j)
                  end if
               end if
               if (av(i+1,j) .eq. 1) then ! western closed boundary
                  select case(Am_method)
                     case (AM_CONSTANT)
                        work2d(i,j)=Am_const*shearX(i,j)
                     case (AM_LES)
                        work2d(i,j)=AmX(i,j)*shearX(i,j)
                  end select
                  if (Am_method.eq.AM_CONSTANT .or. Am_method.eq.AM_LES) then
                     if (calc_phydis) then
                        phydis_wrk(i,j) = work2d(i,j) * dvdxX(i,j)
                     end if
                     work2d(i,j) = DXX * DV(i+1,j) * work2d(i,j)
                  end if
               end if
#endif
#endif
            end if
         end do
#ifdef SLICE_MODEL
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
#endif
         do i=imin,imax          ! VEx defined on V-points
            if (av(i,j).eq.1 .or. av(i,j).eq.2) then
               VEx(i,j)=VEx(i,j)-(work2d(i  ,j)-work2d(i-1,j))*ARVD1
               if (calc_phydis) then
                  phydis_vel(i,j) = _HALF_ * ( phydis_wrk(i-1,j) + phydis_wrk(i,j) )
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
   end if


#ifndef SLICE_MODEL
!  Central for dy(2*Am*dy(V/DV))
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax+1     ! work2d defined on T-points
      do i=imin,imax
         work2d(i,j)=_ZERO_
         if (az(i,j) .eq. 1) then
            select case(Am_method)
               case (AM_LAPLACE)
                  work2d(i,j)=Am_const*dvdyC(i,j)
               case (AM_LES)
                  work2d(i,j)=_TWO_*AmC(i,j)*dvdyC(i,j)
               case (AM_CONSTANT)
                  work2d(i,j)=_TWO_*Am_const*dvdyC(i,j)
            end select
            if (Am_method .ne. NO_AM) then
               if (calc_phydis) then
                  phydis_wrk(i,j) = work2d(i,j) * dvdyC(i,j)
               end if
               work2d(i,j) = DXC * D(i,j) * work2d(i,j)
            end if
            select case(An_method)
               case(1)
                  work2d(i,j)=work2d(i,j)+An_const*DXC*(V(i,j)-V(i,j-1))/DYC
               case(2)
                  if (AnC(i,j) .gt. _ZERO_) then
                     work2d(i,j)=work2d(i,j)+AnC(i,j)*DXC*(V(i,j)-V(i,j-1))/DYC
                  end if
            end select
         else
            if (calc_phydis) then
               phydis_wrk(i,j) = _ZERO_
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax             ! VEx defined on V-points
      do i=imin,imax
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j)=VEx(i,j)-(work2d(i,j+1)-work2d(i,j  ))*ARVD1
            if (calc_phydis) then
               phydis_vel(i,j) = phydis_vel(i,j) + _HALF_*(phydis_wrk(i,j)+phydis_wrk(i,j+1))
            end if
         end if
      end do
   end do
!$OMP END DO
#endif

   if (calc_phydis) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               phydis(i,j) = phydis(i,j) + _HALF_*(phydis_vel(i,j-1)+phydis_vel(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
   end if

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   UEx(imin:imax,j+1) = UEx(imin:imax,j)
   VEx(imin:imax,j-1) = VEx(imin:imax,j)
   VEx(imin:imax,j+1) = VEx(imin:imax,j)
#endif

   if (present(hsd_u)) then
      hsd_u = UEx - hsd_u
   end if
   if (present(hsd_v)) then
      hsd_v = VEx - hsd_v
   end if

#ifdef DEBUG
     write(debug,*) 'Leaving uv_diff_2dh()'
     write(debug,*)
#endif
   return
   end subroutine uv_diff_2dh
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
