#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ww_momentum_3d - continuity eq.\ \label{sec-ww-momentum-3d}
!
! !INTERFACE:
   subroutine ww_momentum_3d()
!
! !DESCRIPTION:
!
! Here, the local continuity equation is calculated in order to obtain
! the grid-related vertical velocity $\bar w_k$. An layer-integrated equation
! for this quantity is given as equation (\ref{ContiLayerInt}) which
! has been derived
! from the differential formulation (\ref{Konti}).
!
! Since the kinematic boundary condition must hold (and is used for the
! derivation of (\ref{ContiLayerInt})), the grid-related vertical
! velocity at the surface muzst be zero, i.e.\ $\bar w_{k_{\max}}=0$.
! This is a good consistence check for the mode splitting, since this is
! only fulfilled if the vertically integrated continuity equation
! (which is the sea surface elevation equation calculated on the
! micro time step) and this local continuity equation are compatible.
!
! The physical vertical velocity is then recalculated from the grid-related
! vertical velocity by means of (\ref{conservative_w}), ... which should
! soon be coded in the routine {\tt tow} in the directory {\tt futils}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,kmin,uu,vv,ww,ho,hn
!  #define CALC_HALO_WW
#ifndef CALC_HALO_WW
   use domain, only: az
   use halo_zones, only: update_3d_halo,wait_halo,z_TAG
#endif
   use getm_timers, only: tic, toc, TIM_WWMOMENTUM, TIM_WWMOMENTUMH
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dtm1
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ww_momentum_3d() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_WWMOMENTUM)

   dtm1=_ONE_/dt

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k)

! OMP-NOTE: k-1 is used at layer k, so we have to conclude
!    one layer at a time (wait after each later).
   do k=1,kmax-1
!$OMP DO SCHEDULE(RUNTIME)
#ifdef CALC_HALO_WW
#ifndef SLICE_MODEL
      do j=jmin-1,jmax+HALO
#endif
         do i=imin-1,imax+HALO
#else
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
#endif
            if (az(i,j) .eq. 1) then
               if (k .lt. kmin(i,j)) then
                  ww(i,j,k)= _ZERO_
               else
                  ww(i,j,k) =   ww(i,j,k-1)                             &
                              - ( hn(i,j,k) - ho(i,j,k) )*dtm1          &
                              - (                                       &
                                   uu(i,j,k)*DYU - uu(i-1,j  ,k)*DYUIM1 &
#ifndef SLICE_MODEL
                                 + vv(i,j,k)*DXV - vv(i  ,j-1,k)*DXVJM1 &
#endif
                                )*ARCD1
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
   end do

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   ww(:,j+1,:) = ww(:,j,:)
#endif

#ifndef CALC_HALO_WW
   call tic(TIM_WWMOMENTUMH)
   call update_3d_halo(ww,ww,az,imin,jmin,imax,jmax,kmax,z_TAG)
   call wait_halo(z_TAG)
   call toc(TIM_WWMOMENTUMH)
#endif

!     Consistency test: ww(i,j,kmax) must always be zero !
   call toc(TIM_WWMOMENTUM)
#ifdef DEBUG
   write(debug,*) 'Leaving ww_momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine ww_momentum_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
