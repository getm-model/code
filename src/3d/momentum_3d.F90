!$Id: momentum_3d.F90,v 1.16 2010-03-25 11:48:55 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: momentum_3d - 3D-momentum for all interior points.
!
! !INTERFACE:
   subroutine momentum_3d(runtype,n)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: kmax,vert_cord,maxdepth
   use variables_3d, only: uu,vv,ho,hn,huo,hun,hvo,hvn
   use m3d, only: ufirst,bdy3d,vert_cord,cord_relax
   use m3d, only: vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,nonhyd_method
   use nonhydrostatic, only: do_nonhydrostatic
   use nonhydrostatic, only: nonhyd_iters
   use variables_3d, only: uu_0,vv_0
   use getm_timers, only: tic,toc,TIM_INTEGR3D,TIM_NH_OVERHEAD
   use internal_pressure, only: do_internal_pressure

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)       :: runtype
   integer,intent(in)       :: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: nonhyd_loop
   logical :: ufirst_0
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'momentum_3d() # ',Ncall
#endif

   if (nonhyd_iters .gt. 1) then
      call tic(TIM_INTEGR3D)
      call tic(TIM_NH_OVERHEAD)
      uu_0  = uu  ; vv_0  = vv
      ufirst_0 = ufirst
      call toc(TIM_NH_OVERHEAD)
      call toc(TIM_INTEGR3D)
   end if
   nonhyd_loop = 1
   do while (nonhyd_loop .le. nonhyd_iters)
      if (nonhyd_loop .gt. 1) then
         call tic(TIM_INTEGR3D)
         call tic(TIM_NH_OVERHEAD)
         uu  = uu_0  ; vv  = vv_0
         ufirst = ufirst_0
         call toc(TIM_INTEGR3D)
      end if
      if (ufirst) then
         call uu_momentum_3d(n,bdy3d)
         call vv_momentum_3d(n,bdy3d)
         ufirst = .false.
      else
         call vv_momentum_3d(n,bdy3d)
         call uu_momentum_3d(n,bdy3d)
         ufirst = .true.
      end if
      if (kmax .gt. 1) then
         call ww_momentum_3d()
      end if
      if (nonhyd_loop .gt. 1) then
         call toc(TIM_NH_OVERHEAD)
      end if
      if (nonhyd_method .ne. 0) then
         call tic(TIM_NH_OVERHEAD)
         call do_nonhydrostatic(nonhyd_loop,vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver)
         call toc(TIM_NH_OVERHEAD)
      end if
      if (nonhyd_loop .lt. nonhyd_iters) then
         call tic(TIM_NH_OVERHEAD)
         call do_internal_pressure()
         call toc(TIM_NH_OVERHEAD)
      end if
      nonhyd_loop = nonhyd_loop + 1
   end do

   if (nonhyd_iters .gt. 1) then
      call tic(TIM_INTEGR3D)
      call tic(TIM_NH_OVERHEAD)
      ufirst = ( .not. ufirst_0 )
      call toc(TIM_NH_OVERHEAD)
      call toc(TIM_INTEGR3D)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine momentum_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
