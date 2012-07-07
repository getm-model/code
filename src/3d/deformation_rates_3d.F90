#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: deformation_rates_3d - \label{deformation_rates_3d}
!
! !INTERFACE:
   subroutine deformation_rates_3d()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: kmax
   use m2d, only: deformation_rates
   use m3d, only: deformCX_3d,deformUV_3d
   use variables_3d, only: uu,vv,hun,hvn
   use variables_3d, only: dudxC_3d,dudxV_3d
#ifndef SLICE_MODEL
   use variables_3d, only: dvdyC_3d,dvdyU_3d
#endif
   use variables_3d, only: shearX_3d,shearU_3d
   use getm_timers, only: tic,toc,TIM_DEFORM3D

   IMPLICIT NONE

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
   write(debug,*) 'deformation_rates_3d() # ',Ncall
#endif
   call tic(TIM_DEFORM3D)

   if (deformCX_3d) then
      if (deformUV_3d) then
         do k=1,kmax
            call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k),     &
                                   dudxC=dudxC_3d(:,:,k),dudxV=dudxV_3d(:,:,k),   &
#ifndef SLICE_MODEL
                                   dvdyC=dvdyC_3d(:,:,k),dvdyU=dvdyU_3d(:,:,k),   &
#endif
                                   shearX=shearX_3d(:,:,k),shearU=shearU_3d(:,:,k))
         end do
      else
         do k=1,kmax
            call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k), &
                                   dudxC=dudxC_3d(:,:,k),                     &
#ifndef SLICE_MODEL
                                   dvdyC=dvdyC_3d(:,:,k),                     &
#endif
                                   shearX=shearX_3d(:,:,k))
         end do
      end if
   end if

   call toc(TIM_DEFORM3D)
#ifdef DEBUG
   write(debug,*) 'Leaving deformation_rates_3d()'
   write(debug,*)
#endif
   return
   end subroutine deformation_rates_3d

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
