!$Id: calc_uvex.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: calc_uvex
!
! !INTERFACE:
   subroutine calc_uvex(U,V,D,DU,DV)
!
! !DESCRIPTION:
!
! !USES:
   use m2d, only: Am_method,An_method,NO_AM,AM_CONSTANT,AM_LES
   use variables_2d, only: UEx,VEx
   use variables_2d, only: dudxC,dudxV,shearX,shearU
#ifndef SLICE_MODEL
   use variables_2d, only: dvdyC,dvdyU
#endif
   use les, only: do_les_2d
   use variables_les, only: AmC_2d,AmX_2d
   use getm_timers,  only: tic,toc,TIM_UVEX
!  needed by interface headers!
   use domain, only: imin,imax,jmin,jmax

   IMPLICIT NONE

#include "deformation_rates.h"
#include "uv_diffusion.h"

!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: U,V,D,DU,DV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'calc_uvex() # ',Ncall
#endif
   CALL tic(TIM_UVEX)

      call uv_advect(U,V,DU,DV)
      select case(Am_method)
         case(NO_AM)
            if (An_method .gt. 0) then
               call uv_diffusion(UEx,VEx,U=U,V=V)
            end if
         case(AM_CONSTANT)
            call deformation_rates(U,V,DU,DV,   &
                                   dudxC=dudxC, &
#ifndef SLICE_MODEL
                                   dvdyC=dvdyC, &
#endif
                                   shearX=shearX)
            call uv_diffusion(UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV, &
                              dudxC=dudxC,                     &
#ifndef SLICE_MODEL
                              dvdyC=dvdyC,                     &
#endif
                              shearX=shearX)
         case(AM_LES)
            call deformation_rates(U,V,DU,DV,                 &
                                   dudxC=dudxC,dudxV=dudxV,   &
#ifndef SLICE_MODEL
                                   dvdyC=dvdyC,dvdyU=dvdyU,   &
#endif
                                   shearX=shearX,shearU=shearU)
            CALL toc(TIM_UVEX)
            call do_les_2d(dudxC,dudxV, &
#ifndef SLICE_MODEL
                           dvdyC,dvdyU, &
#endif
                           shearX,shearU)
            CALL tic(TIM_UVEX)
            call uv_diffusion(UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV,   &
                              dudxC=dudxC,                       &
#ifndef SLICE_MODEL
                              dvdyC=dvdyC,                       &
#endif
                              shearX=shearX,AmC=AmC_2d,AmX=AmX_2d)
      end select

   CALL toc(TIM_UVEX)
#ifdef DEBUG
   write(debug,*) 'Leaving calc_uvex()'
   write(debug,*)
#endif
   return
   end subroutine calc_uvex

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
