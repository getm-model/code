#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: calc_uvex
!
! !INTERFACE:
   subroutine calc_uvex(An_method,U,V,D,DU,DV)
!  Note (KK): keep in sync with interface in m2d_general.F90
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use m2d_general, only: uv_diffusion
   use m2d, only: Am
   use variables_2d, only: UEx,VEx
   use getm_timers,  only: tic,toc,TIM_UVEX

   IMPLICIT NONE

!
! !INPUT PARAMETERS:
   integer,intent(in)                      :: An_method
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

   UEx=_ZERO_ ; VEx=_ZERO_

#ifdef NO_ADVECT
   STDERR 'NO_ADVECT 2D'
#else
#ifndef UV_ADV_DIRECT
   call uv_advect(U,V,DU,DV)
   if (Am.gt._ZERO_ .or. An_method.gt.0) then
      call uv_diffusion(An_method,UEx,VEx,U=U,V=V,D=D,DU=DU,DV=DV)
   end if
#endif
#endif

   CALL toc(TIM_UVEX)
#ifdef DEBUG
   write(debug,*) 'Leaving calc_uvex()'
   write(debug,*)
#endif
   return
   end subroutine calc_uvex

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
