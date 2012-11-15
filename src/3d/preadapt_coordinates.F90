#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  preadapt_coordinates - pre-adaptation of the vertical coordinate
! \label{sec-preadapt-coordinates}
!
! !INTERFACE:
   subroutine preadapt_coordinates(preadapt)
!
! !DESCRIPTION:
!
! Pre-adaptation of the vertical coordinates using stratification
! and water depth. The initial salinities and temperatures are re-interpolated
! onto the new grid after pre-adaptation
!
! !USES:
   use getm_timers, only: tic, toc, TIM_COORDS
#ifndef NO_BAROCLINIC
   use m3d, only: calc_salt, calc_temp
   use salinity, only: init_salinity_field,
   use temperature, only: init_temperature_field
   use eqstate, only: do_eqstate
#endif
   use variables_3d, only: SS
   use domain, only: vert_cord
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: preadapt
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
    integer         :: ii
!EOP
!-----------------------------------------------------------------------
!BOC

   call tic(TIM_COORDS)

   if (preadapt.ne.0) then

      LEVEL1 'pre-adapting coordinates'

      do ii=1,preadapt

#ifndef NO_BAROCLINIC
         call do_eqstate()
         call buoyancy_frequency()
#endif
         call adaptive_coordinates(.false.,.false.)
#ifndef NO_BAROCLINIC
         if(calc_salt) then
            LEVEL2 'reinterpolating initial salinity'
            call init_salinity_field()
         end if
         if(calc_temp) then
            LEVEL2 'reinterpolating initial temperature'
            call init_temperature_field()
         end if
         call do_eqstate()
         call buoyancy_frequency()
#endif

      end do

   end if

   call toc(TIM_COORDS)
#ifdef DEBUG
   write(debug,*) 'Leaving preadapt_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine preadapt_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------

