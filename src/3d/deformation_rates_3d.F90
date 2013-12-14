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
   use domain, only: imin,imax,jmin,jmax,kmax
   use pool, only: deformation_rates
   use variables_3d, only: deformC_3d,deformX_3d,deformUV_3d
   use variables_3d, only: uu,vv,hun,hvn
   use variables_3d, only: dudxC_3d,dudxV_3d,dvdxX_3d
   use variables_3d, only: dvdyC_3d,dvdyU_3d,dudyX_3d
   use variables_3d, only: shearX_3d,shearU_3d
   use variables_3d, only: do_numerical_analyses_3d
   use getm_timers, only: tic,toc,TIM_DEFORM3D

   IMPLICIT NONE

!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   type t_pa2d
      REALTYPE,dimension(:,:),pointer :: p2d
   end type t_pa2d
   type(t_pa2d),dimension(1:kmax)      :: pa_dvdxX,pa_dudyX
   logical                             :: calc_dvdxX,calc_dudyX
#ifndef _POINTER_REMAP_
   REALTYPE,dimension(I2DFIELD),target :: dvdxX,dudyX
#endif
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

   calc_dvdxX = associated(dvdxX_3d)
   if (calc_dvdxX) then
      do k=1,kmax
#ifdef _POINTER_REMAP_
         pa_dvdxX(k)%p2d(imin-HALO:,jmin-HALO:) => dvdxX_3d(:,:,k)
#else
         pa_dvdxX(k)%p2d => dvdxX
#endif
      end do
   else
      do k=1,kmax
         pa_dvdxX(k)%p2d => null()
      end do
   end if

   calc_dudyX = associated(dudyX_3d)
   if (calc_dudyX) then
      do k=1,kmax
#ifdef _POINTER_REMAP_
         pa_dudyX(k)%p2d(imin-HALO:,jmin-HALO:) => dudyX_3d(:,:,k)
#else
         pa_dudyX(k)%p2d => dudyX
#endif
      end do
   else
      do k=1,kmax
         pa_dudyX(k)%p2d => null()
      end do
   end if

   if (deformC_3d) then
      if (deformX_3d) then
         if (do_numerical_analyses_3d) then
            if (deformUV_3d) then
               do k=1,kmax
                  call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k),                         &
                                         dudxC=dudxC_3d(:,:,k),dudxV=dudxV_3d(:,:,k),dvdxX=pa_dvdxX(k)%p2d, &
#ifndef SLICE_MODEL
                                         dvdyC=dvdyC_3d(:,:,k),dvdyU=dvdyU_3d(:,:,k),dudyX=pa_dudyX(k)%p2d, &
#endif
                                         shearX=shearX_3d(:,:,k),shearU=shearU_3d(:,:,k))
#ifndef _POINTER_REMAP_
                  if (calc_dvdxX) then
                     dvdxX_3d(:,:,k) = pa_dvdxX(k)%p2d
                  end if
                  if (calc_dudyX) then
                     dudyX_3d(:,:,k) = pa_dudyX(k)%p2d
                  end if
#endif
               end do
            else
               do k=1,kmax
                  call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k),   &
                                         dudxC=dudxC_3d(:,:,k),dvdxX=pa_dvdxX(k)%p2d, &
#ifndef SLICE_MODEL
                                         dvdyC=dvdyC_3d(:,:,k),dudyX=pa_dudyX(k)%p2d, &
#endif
                                         shearX=shearX_3d(:,:,k))
#ifndef _POINTER_REMAP_
                  if (calc_dvdxX) then
                     dvdxX_3d(:,:,k) = pa_dvdxX(k)%p2d
                  end if
                  if (calc_dudyX) then
                     dudyX_3d(:,:,k) = pa_dudyX(k)%p2d
                  end if
#endif
               end do
            end if
         else
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
      else
         do k=1,kmax
            call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k), &
                                   dudxC=dudxC_3d(:,:,k)                      &
#ifndef SLICE_MODEL
                                  ,dvdyC=dvdyC_3d(:,:,k)                      &
#endif
                                  )
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
