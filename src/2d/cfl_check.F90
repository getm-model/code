#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cfl_check - check for explicit barotropic time step.
!
! !INTERFACE:
   subroutine cfl_check()
!
! !DESCRIPTION:
!
! This routine loops over all horizontal grid points and calculated the
! maximum time step according to the shallow water criterium by
! \cite{BECKERSea93}:
!
! \begin{equation}
! \Delta t_{\max} = \min_{i,j} \left\{\frac{\Delta x_{i,j} \Delta y_{i,j}}
! {\sqrt{2} c_{i,j} \sqrt{\Delta x_{i,j}^2+ \Delta y_{i,j}^2}}\right\}
! \end{equation}
!
! with the local Courant number
!
! \begin{equation}
! c_{i,j}=\sqrt{g H_{i,j}},
! \end{equation}
!
! where $g$ is the gravitational acceleration and $H_{i,j}$ is the local
! bathymetry value. In case that the chosen micro time step $\Delta t_m$
! is larger than $\Delta t_{\max}$, the program will be aborted. In any
! the CFL diagnostics will be written to standard output.
!

!
! !USES:
   use parameters, only: g
   use domain, only: imin,imax,jmin,jmax,H,az
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,dxc
#else
   use domain, only: dy,dx
#endif
   use variables_2d, only: dtm
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: pos(2),max_pos(2),rc,i,j
   REALTYPE                  :: h_max=-99.,c,max_dt,dtt,dxeff
   logical, dimension(:,:), allocatable :: lmask
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'CFL check'

   allocate(lmask(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'cfl_check: Error allocating memory (lmask)'

   lmask = .false.
   lmask = (az .gt. 0)
   h_max = maxval(H,mask = lmask)
   pos = maxloc(H,mask = lmask)
#if defined(SPHERICAL) || defined(CURVILINEAR)
   max_dt=1000000000.
   do i=imin,imax
      do j=jmin,jmax
         if (az(i,j) .ge. 1 .and. H(i,j) .gt. _ZERO_) then
#if 0
            dtt=min(dxc(i,j),dyc(i,j))/sqrt(2.*g*H(i,j))
#else
            c = sqrt(g*H(i,j))
#ifdef SLICE_MODEL
            dxeff = dxc(i,j)
#else
            dxeff = (dxc(i,j)*dyc(i,j))/ &
                     (sqrt(2.0)*sqrt(dxc(i,j)*dxc(i,j)+dyc(i,j)*dyc(i,j)))
#endif
            dtt = dxeff/c

#endif
            if (dtt .lt. max_dt) then
               max_dt=dtt
               max_pos(1)=i
               max_pos(2)=j
            end if
         end if
      end do
   end do
   LEVEL3 'max CFL number at depth=',H(max_pos(1),max_pos(2)),' at ',max_pos
   LEVEL3 'at this position, dx = ',dxc(max_pos(1),max_pos(2)),' and dy =  ',dyc(max_pos(1),max_pos(2))
#else
   c = sqrt(g*h_max)
!  Becker and Deleersnijder
#ifdef SLICE_MODEL
   dxeff = dx
#else
   dxeff = (dx*dy)/(sqrt(2.0)*sqrt(dx*dx+dy*dy))
#endif
   max_dt = dxeff/c

#if 0
#ifdef POM_CFL
   max_dt=0.5/(c*sqrt( _ONE_/dx**2 + _ONE_/dy**2 ))
#endif
#ifdef OLD_GETM
   max_dt = min(dx,dy)/(c*sqrt(2.))
#endif
#endif

   LEVEL3 'max depth=',h_max,' at ',pos
#endif

   if (dtm .gt. max_dt) then
      FATAL 'reduce time-step (',dtm,') below ',max_dt
      stop 'cfl_check()'
   else
      LEVEL3 'used dt (',dtm,') is less than ',max_dt
   end if

#ifdef FORTRAN90
   deallocate(lmask,stat=rc)
   if (rc /= 0) stop 'cfl_check: Error allocating memory (lmask)'
#endif
   return
   end subroutine cfl_check
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
