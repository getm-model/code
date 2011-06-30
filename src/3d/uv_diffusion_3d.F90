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
   use m2d, only: Am_method
   use variables_3d, only: uuEx,vvEx,hn,hun,hvn
   use variables_3d, only: dudxC_3d,dvdyc_3d,shearX_3d
#ifdef _LES_
   use variables_les, only: Am_3d,AmX_3d
#endif
   use getm_timers, only: tic, toc, TIM_UVDIFF3D

   IMPLICIT NONE

#include "../2d/uv_diffusion.h"
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: k

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion_3d() # ',Ncall
#endif
   call tic(TIM_UVDIFF3D)

!  KK-TODO: correct timers
   do k=1,kmax
      select case(Am_method)
         case(1)
            call uv_diffusion(1,0,uuEx(:,:,k),vvEx(:,:,k),                 &
                              dudxC=dudxC_3d(:,:,k),dvdyC=dvdyC_3d(:,:,k), &
                              shearX=shearX_3d(:,:,k),                     &
                              D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k))
#ifdef _LES_
         case(2)
            call uv_diffusion(2,0,uuEx(:,:,k),vvEx(:,:,k),                 &
                              dudxC=dudxC_3d(:,:,k),dvdyC=dvdyC_3d(:,:,k), &
                              shearX=shearX_3d(:,:,k),                     &
                              D=hn(:,:,k),DU=hun(:,:,k),DV=hvn(:,:,k),     &
                              Am=Am_3d(:,:,k),AmX=AmX_3d(:,:,k))
#endif
      end select
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
