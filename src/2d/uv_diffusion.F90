#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion - lateral diffusion of depth-averaged velocity
! \label{sec-uv-diffusion}
!
! !INTERFACE:
   subroutine uv_diffusion(An_method,U,V,D,DU,DV)

!  Note (KK): keep in sync with interface in m2d.F90
!
! !DESCRIPTION:
! This wrapper calls routine {\tt uv\_diff\_2dh} (see section
! \ref{sec-uv-diff-2dh} on page \pageref{sec-uv-diff-2dh}).
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use m2d, only: deformation_rates,uv_diff_2dh
   use m2d, only: Am_method,NO_AM,AM_LAPLACE,AM_LES,AM_CONSTANT
   use variables_2d, only: UEx,VEx
   use variables_2d, only: dudxC,dudxV,dvdyC,dvdyU,shearX,shearU
   use les, only: do_les_2d
   use variables_les, only: AmC_2d,AmX_2d
   use getm_timers,  only: tic,toc,TIM_UVDIFF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                             :: An_method
   REALTYPE,dimension(E2DFIELD),intent(in)        :: U,V,D
   REALTYPE,dimension(E2DFIELD),target,intent(in) :: DU,DV
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_diffusion() # ',Ncall
#endif
   CALL tic(TIM_UVDIFF)

   select case(Am_method)
      case(NO_AM)
         if (An_method .gt. 0) then
            call uv_diff_2dh(An_method,UEx,VEx,U=U,V=V)
         end if
      case(AM_LAPLACE)
         call deformation_rates(U,V,DU,DV,  &
                                dudxC=dudxC &
#ifndef SLICE_MODEL
                               ,dvdyC=dvdyC &
#endif
                               )
         call uv_diff_2dh(An_method,UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV, &
                          dudxC=dudxC                                &
#ifndef SLICE_MODEL
                         ,dvdyC=dvdyC                                &
#endif
                         )
      case(AM_LES)
         call deformation_rates(U,V,DU,DV,                 &
                                dudxC=dudxC,dudxV=dudxV,   &
#ifndef SLICE_MODEL
                                dvdyC=dvdyC,dvdyU=dvdyU,   &
#endif
                                shearX=shearX,shearU=shearU)
         CALL toc(TIM_UVDIFF)
         call do_les_2d(dudxC,dudxV, &
#ifndef SLICE_MODEL
                        dvdyC,dvdyU, &
#endif
                        shearX,shearU)
         CALL tic(TIM_UVDIFF)
         call uv_diff_2dh(An_method,UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV, &
                          dudxC=dudxC,                               &
#ifndef SLICE_MODEL
                          dvdyC=dvdyC,                               &
#endif
                          shearX=shearX,AmC=AmC_2d,AmX=AmX_2d)
      case(AM_CONSTANT)
         call deformation_rates(U,V,DU,DV,   &
                                dudxC=dudxC, &
#ifndef SLICE_MODEL
                                dvdyC=dvdyC, &
#endif
                                shearX=shearX)
         call uv_diff_2dh(An_method,UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV, &
                          dudxC=dudxC,                               &
#ifndef SLICE_MODEL
                          dvdyC=dvdyC,                               &
#endif
                          shearX=shearX)

   end select

   CALL toc(TIM_UVDIFF)
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
