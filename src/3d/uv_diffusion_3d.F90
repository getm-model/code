#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion_3d - lateral diffusion of 3D velocity
! \label{sec-uv-diffusion-3d}
!
! !INTERFACE:
   subroutine uv_diffusion_3d()
!
! !DESCRIPTION:
! This wrapper calls routine {\tt uv\_diff\_2dh} (see section
! \ref{sec-uv-diff-2dh} on page \pageref{sec-uv-diff-2dh}) for each
! layer.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use m2d, only: uv_diff_2dh
   use m2d, only: Am_method,AM_LAPLACE,AM_LES,AM_CONSTANT
   use variables_3d, only: uu,vv,uuEx,vvEx,hn,hun,hvn
   use variables_3d, only: dudxC_3d,dvdyC_3d
   use variables_3d, only: dvdxX_3d,dudyX_3d,shearX_3d
   use variables_3d, only: do_numerical_analyses_3d
   use variables_3d, only: phydis_3d
   use variables_les, only: AmC_2d,AmX_2d,AmC_3d,AmX_3d
#ifdef _MOMENTUM_TERMS_
   use domain, only: dry_u,dry_v
   use variables_3d, only: hsd_u,hsd_v
#endif
   use getm_timers, only: tic, toc, TIM_UVDIFF3D
!$ use omp_lib
   IMPLICIT NONE

!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   type t_pa2d
      REALTYPE,dimension(:,:),pointer :: p2d
   end type t_pa2d
   type(t_pa2d),dimension(1:kmax)      :: pa_pd2d,pa_dvdxX,pa_dudyX
   logical                             :: calc_phydis
#ifndef _POINTER_REMAP_
   REALTYPE,dimension(I2DFIELD),target :: pd2d,dvdxX,dudyX
#endif
   integer :: i,j,k

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion_3d() # ',Ncall
#endif
   call tic(TIM_UVDIFF3D)

   calc_phydis = associated(phydis_3d)
   if (calc_phydis) then
      do k=1,kmax
#ifdef _POINTER_REMAP_
         pa_pd2d (k)%p2d(imin-HALO:,jmin-HALO:) => phydis_3d(:,:,k)
         pa_dvdxX(k)%p2d(imin-HALO:,jmin-HALO:) => dvdxX_3d (:,:,k)
#ifndef SLICE_MODEL
         pa_dudyX(k)%p2d(imin-HALO:,jmin-HALO:) => dudyX_3d (:,:,k)
#endif
#else
         pa_pd2d (k)%p2d => pd2d
         pa_dvdxX(k)%p2d => dvdxX
         pa_dudyX(k)%p2d => dudyX
#endif
      end do
   else
      do k=1,kmax
         pa_pd2d(k)%p2d => null()
      end do
   end if

   select case(Am_method)
      case(AM_LAPLACE)
         do k=1,kmax
            call uv_diff_2dh(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                             D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),           &
                             dudxC=dudxC_3d(:,:,k),                             &
#ifndef SLICE_MODEL
                             dvdyC=dvdyC_3d(:,:,k),                             &
#endif
                             phydis=pa_pd2d(k)%p2d                              &
#ifdef _MOMENTUM_TERMS_
                            ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)              &
#endif
                            )
#ifndef _POINTER_REMAP_
            if (calc_phydis) then
               phydis_3d(:,:,k) = pa_pd2d(k)%p2d
            end if
#endif
         end do
      case(AM_LES)
         do k=1,kmax
#ifndef _POINTER_REMAP_
            if (calc_phydis) then
               pa_dvdxX(k)%p2d = dvdxX_3d(:,:,k)
#ifndef SLICE_MODEL
               pa_dudyX(k)%p2d = dudyX_3d(:,:,k)
#endif
            end if
#endif
            call uv_diff_2dh(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                             D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),           &
                             dudxC=dudxC_3d(:,:,k),dvdxX=pa_dvdxX(k)%p2d,       &
#ifndef SLICE_MODEL
                             dvdyC=dvdyC_3d(:,:,k),dudyX=pa_dudyX(k)%p2d,       &
#endif
                             shearX=shearX_3d(:,:,k),                           &
                             AmC=AmC_3d(:,:,k),AmX=AmX_3d(:,:,k),               &
                             phydis=pa_pd2d(k)%p2d                              &
#ifdef _MOMENTUM_TERMS_
                            ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)              &
#endif
                            )
#ifndef _POINTER_REMAP_
            if (calc_phydis) then
               phydis_3d(:,:,k) = pa_pd2d(k)%p2d
            end if
#endif
         end do
      case(AM_CONSTANT)
         do k=1,kmax
#ifndef _POINTER_REMAP_
            if (calc_phydis) then
               pa_dvdxX(k)%p2d = dvdxX_3d(:,:,k)
#ifndef SLICE_MODEL
               pa_dudyX(k)%p2d = dudyX_3d(:,:,k)
#endif
            end if
#endif
            call uv_diff_2dh(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                             D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),           &
                             dudxC=dudxC_3d(:,:,k),dvdxX=pa_dvdxX(k)%p2d,       &
#ifndef SLICE_MODEL
                             dvdyC=dvdyC_3d(:,:,k),dudyX=pa_dudyX(k)%p2d,       &
#endif
                             shearX=shearX_3d(:,:,k),                           &
                             phydis=pa_pd2d(k)%p2d                              &
#ifdef _MOMENTUM_TERMS_
                            ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)              &
#endif
                            )
#ifndef _POINTER_REMAP_
            if (calc_phydis) then
               phydis_3d(:,:,k) = pa_pd2d(k)%p2d
            end if
#endif
         end do
   end select

#ifdef _MOMENTUM_TERMS_
!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k)

   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            hsd_u(i,j,k) = dry_u(i,j) * hsd_u(i,j,k)
            hsd_v(i,j,k) = dry_v(i,j) * hsd_v(i,j,k)
         end do
      end do
!$OMP END DO NOWAIT
   end do

!$OMP END PARALLEL
#endif

   call toc(TIM_UVDIFF3D)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_diffusion_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_diffusion_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
