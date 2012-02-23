#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion_3d - hor.\ momentum diffusion
! \label{sec-uv-diffusion-3d}
!
! !INTERFACE:
   subroutine uv_diffusion_3d()
!
! !DESCRIPTION:
!  This wrapper calls uv_diffusion for each layer.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use m2d_general, only: uv_diffusion
   use variables_3d, only: uu,vv,uuEx,vvEx,hn,hun,hvn
#ifdef _MOMENTUM_TERMS_
   use domain, only: dry_u,dry_v
   use variables_3d, only: hsd_u,hsd_v
#endif
   use getm_timers, only: tic, toc, TIM_UVDIFF3D
!$ use omp_lib
   IMPLICIT NONE

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

   do k=1,kmax
      call uv_diffusion(0,uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                        D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k)            &
#ifdef _MOMENTUM_TERMS_
                        ,hsd_u=hsd_u(:,:,k),hsd_v=hsd_v(:,:,k)             &
#endif
                       )
   end do

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
