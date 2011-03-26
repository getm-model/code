#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stop_macro - terminates the macro loop \label{sec-stop-macro}
!
! !INTERFACE:
   subroutine stop_macro
!
! !DESCRIPTION:
!
! This routine should be called from {\tt m3d} at the end of each macro
! time step in order to copy the vertically interated and temporally
! averaged transports to old values {\tt Uinto} and {\tt Vinto}, and
! to reinitialise the transports {\tt Uint} and {\tt Vint}
! to zero.
! This is also the place where the surface velocity divergence is
! calculated, i.e.\ {\tt div}$=\partial_x u(z=\eta)+\partial_y v(z=\eta)$,
! to be output in the 2d netcdf file.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use variables_2d, only: Uint,Uinto,Vint,Vinto,surfdiv
   use variables_3d, only: hun,hvn,uu,vv
   use getm_timers, only: tic, toc, TIM_STOPMCR
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
   use domain, only: dx,dy,ard1
#endif
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stop_macro() # ',Ncall
#endif

  Uinto=Uint
  Uint= _ZERO_
  Vinto=Vint
  Vint= _ZERO_

#ifdef DEBUG
   write(debug,*) 'Leaving stop_macro()'
   write(debug,*)
#endif
   return
   end subroutine stop_macro
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
