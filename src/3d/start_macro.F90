!$Id: start_macro.F90,v 1.12 2009-08-18 10:24:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: start_macro - initialise the macro loop \label{sec-start-macro}
!
! !INTERFACE:
   subroutine start_macro()
!
! !DESCRIPTION:
!
! This routine needs to be called from {\tt m3d} at the beginning
! of each macro time step. Here, the sea surface elevations at the
! before and after the macro time step are updated at the
! T-, U- and V-points.the sea surface elevations at the
! before and after the macro time step are updated at the
! T-, U- and V-points, their notation is {\tt sseo}, {\tt ssen},
! {\tt ssuo}, {\tt ssun}, {\tt ssvo} and {\tt ssvn}, where {\tt e},
! {\tt u} and {\tt v} stand for T-, U- and V-point and {\tt o} and
! {\tt n} for old and new, respectively, see also the description of
! {\tt variables\_3d} in section \ref{sec-variables-3d} on page 
! \pageref{sec-variables-3d}.
!
! Furthermore, the vertically integrated transports {\tt Uint}
! and {\tt Vint} are here divided by the number of micro time
! steps per macro time step, {\tt M}, in order to obtain 
! the time-averaged transports.
!
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,H,HU,HV,min_depth
   use m2d, only: z,Uint,Vint
   use m3d, only: M
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   use getm_timers, only: tic, toc, TIM_STARTMCR
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: split
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'start_macro() # ',Ncall
#endif
   call tic(TIM_STARTMCR)

   do j=jmin-HALO,jmax+HALO         ! Defining 'old' and 'new' sea surface
      do i=imin-HALO,imax+HALO      ! elevation for macro time step
         sseo(i,j)=ssen(i,j)
         ssen(i,j)=z(i,j)
      end do
   end do

   do j=jmin-HALO,jmax+HALO             ! Same for U-points
      do i=imin-HALO,imax+HALO-1
         ssuo(i,j)=ssun(i,j)
         ssun(i,j)=0.25*(sseo(i,j)+sseo(i+1,j)+ssen(i,j)+ssen(i+1,j))
         ssun(i,j)=max(ssun(i,j),-HU(i,j)+min_depth)
      end do
   end do

   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO             ! Same for V-points
         ssvo(i,j)=ssvn(i,j)
         ssvn(i,j)=0.25*(sseo(i,j)+sseo(i,j+1)+ssen(i,j)+ssen(i,j+1))
         ssvn(i,j)=max(ssvn(i,j),-HV(i,j)+min_depth)
      end do
   end do

! Defining vertically integrated, conservative
! u- and v-transport for macro time step

   split = _ONE_/float(M)
   Uint = split*Uint
   Vint = split*Vint

   call toc(TIM_STARTMCR)
#ifdef DEBUG
   write(debug,*) 'Leaving start_macro()'
   write(debug,*)
#endif
   return
   end subroutine start_macro
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
