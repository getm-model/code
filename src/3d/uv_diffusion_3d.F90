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
   use m2d, only: Am_method,AM_CONSTANT,AM_LES
   use variables_3d, only: uu,vv,uuEx,vvEx,hn,hun,hvn
   use variables_3d, only: dudxC_3d,dvdyC_3d,shearX_3d
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

   select case(Am_method)
      case(AM_CONSTANT)
         do k=1,kmax
            call uv_diff_2dh(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                             D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),           &
                             dudxC=dudxC_3d(:,:,k),                             &
#ifndef SLICE_MODEL
                             dvdyC=dvdyC_3d(:,:,k),                             &
#endif
                             shearX=shearX_3d(:,:,k)                            &
#ifdef _MOMENTUM_TERMS_
                             ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)             &
#endif
                            )
         end do
      case(AM_LES)
         do k=1,kmax
            call uv_diff_2dh(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                             D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),           &
                             dudxC=dudxC_3d(:,:,k),                             &
#ifndef SLICE_MODEL
                             dvdyC=dvdyC_3d(:,:,k),                             &
#endif
                             shearX=shearX_3d(:,:,k),                           &
                             AmC=AmC_3d(:,:,k),AmX=AmX_3d(:,:,k)                &
#ifdef _MOMENTUM_TERMS_
                             ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)             &
#endif
                            )
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
