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
!
!
! !USES:
   use domain, only: kmax
   use m2d, only: Am_method,AM_CONSTANT,AM_LES
   use variables_3d, only: uu,vv,uuEx,vvEx,hn,hun,hvn
   use variables_3d, only: dudxC_3d,shearX_3d
#ifndef SLICE_MODEL
   use variables_3d, only: dvdyc_3d
#endif
   use variables_les, only: AmC_2d,AmX_2d,AmC_3d,AmX_3d
   use getm_timers, only: tic, toc, TIM_UVDIFF3D

   IMPLICIT NONE

#include "../2d/uv_diffusion.h"
!
! !INPUT PARAMETERS:
!
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

   select case(Am_method)
      case(AM_CONSTANT)
         do k=1,kmax
            call uv_diffusion(uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                              D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),         &
                              dudxC=dudxC_3d(:,:,k),                           &
#ifndef SLICE_MODEL
                              dvdyC=dvdyC_3d(:,:,k),                           &
#endif
                              shearX=shearX_3d(:,:,k))
         end do
      case(AM_LES)
         do k=1,kmax
            call uv_diffusion(uuEx(:,:,k),vvEx(:,:,k),U=uu(:,:,k),V=vv(:,:,k), &
                              D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),         &
                              dudxC=dudxC_3d(:,:,k),                           &
#ifndef SLICE_MODEL
                              dvdyC=dvdyC_3d(:,:,k),                           &
#endif
                              shearX=shearX_3d(:,:,k),                         &
                              AmC=AmC_3d(:,:,k),AmX=AmX_3d(:,:,k))
         end do
   end select

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
