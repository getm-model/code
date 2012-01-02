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
   use domain, only: kmax
   use m2d_general, only: uv_diffusion
   use variables_3d, only: uu,vv,uuEx,vvEx,hn,hun,hvn
   use getm_timers, only: tic, toc, TIM_UVDIFF3D

   IMPLICIT NONE

! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: k

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
                        D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k))
   end do

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
