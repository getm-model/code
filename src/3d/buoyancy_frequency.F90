#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: buoyancy_frequency - calculates the buoyancy frequency
! \label{sec-buoyancy-frequency}
!
! !INTERFACE:
   subroutine buoyancy_frequency()
!
! !DESCRIPTION:
!
! Here, the buoyancy frequency squared, $N^2=\partial_z b$, with buoyancy $b$ from
! (\ref{bdef}) is calculated.
! Two alternative methods are coded.
! The straight-forward method which is explained first,
! has the disadvantage of generating numerical instabilities.
!
! The straight-forward discretisation of
! $N^2$ is given by
!
! \begin{equation}\label{Nstraight}
! \begin{array}{l}
! \displaystyle
! (N^2)_{i,j,k}\approx
! \frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}.
! \end{array}
! \end{equation}
!
! In some cases, together with the straight-forward discretisation
! of the shear squared, (\ref{ShearSquaredOld}), this
! did not produce stable numerical results. The reason for this might be that
! the velocities involved in the calculation for the shear squared do depend
! on the buoyancies in the two
! neighbouring T-points such that the straight-forward
! method (\ref{Nstraight}) leads to an inconsistency.
! However, other experiments with the energy-conserving discretisation
! of the shear stress squared, (\ref{ShearSquaredNew})
! and the straight-forward discretisation of
! $N^2$, (\ref{Nstraight}),  produced numerically stable results.
!
! Most stable results have been obtained with a weighted average
! for the $N^2$ calculation:
!
! \begin{equation}\label{Naveraged}
! \begin{array}{l}
! \displaystyle
! (N^2)_{i,j,k}\approx
! \frac16 \Bigg(
! 2\frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}
! \\ \\ \qquad\qquad\qquad
! \displaystyle
! +
! \frac{b_{i+1,j,k+1}-b_{i+1,j,k}}{\frac12(h^t_{i+1,j,k+1}+h^t_{i+1,j,k})}
! +
! \frac{b_{i-1,j,k+1}-b_{i-1,j,k}}{\frac12(h^t_{i-1,j,k+1}+h^t_{i-1,j,k})}
! \\ \\ \qquad\qquad\qquad
! \displaystyle
! +
! \frac{b_{i,j+1,k+1}-b_{i,j+1,k}}{\frac12(h^t_{i,j+1,k+1}+h^t_{i,j+1,k})}
! +
! \frac{b_{i,j-1,k+1}-b_{i,j-1,k}}{\frac12(h^t_{i,j-1,k+1}+h^t_{i,j-1,k})}
! \Bigg).
! \end{array}
! \end{equation}
!
! These stability issues need to be further investigated in the future.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,au,av,az
   use variables_3d, only: hn
   use parameters, only: g,rho_0
   use variables_3d, only: NN,buoy,T,S
#ifndef _OLD_BVF_
   use variables_3d, only: alpha,beta
#endif
   use getm_timers, only: tic, toc, TIM_NN
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,nb
   REALTYPE                  :: dz,NNc,ttt
   REALTYPE                  :: NNe,NNw,NNn,NNs
   REALTYPE, parameter       :: small_bvf = 1.d-10
#ifdef _SMOOTH_BVF_VERT_
   REALTYPE                  :: below,center,above
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'buoyancy_frequency() # ',Ncall
#endif
   call tic(TIM_NN)

! BJB-TODO: There are a lot of single- (or lower-) precision constants
!   in this routine. Enough that it significantly influences the result
!   if they are correctly converted to full (double) precision.
!   A conversion should be done, but I will not do it as part of
!   the OMP-implementation as it muddles the comparisons of old and
!   new results. BJB 2009-09-22

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,nb,dz,NNc,NNe,NNw,NNn,NNs)

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax+1
      do i=imin-1,imax+1
         if (az(i,j) .ge. 1 ) then
            NN(i,j,kmax) = small_bvf
            do k=kmax-1,1,-1
               dz=_HALF_*(hn(i,j,k+1)+hn(i,j,k))
#ifdef _OLD_BVF_
               NNc =(buoy(i,j,k+1)-buoy(i,j,k))/dz
#else
               NNc = -g / rho_0
               NNc = NNc/dz *(alpha(i,j,k) *(T(i,j,k+1)-T(i,j,k)) &
                         + beta(i,j,k) *(S(i,j,k+1)-S(i,j,k)))
#endif
               if (abs(NNc) .lt. small_bvf ) then
                  NNc = sign(_ONE_,NNc) * small_bvf
               end if
               NN(i,j,k)= NNc
            end do
#ifdef _SMOOTH_BVF_VERT_
            if ( kmax .ge. 4 ) then
               below  = NN(i,j,1)
               center = NN(i,j,2)
               above  = NN(i,j,3)
               do k= 2,kmax-2
                  center = _HALF_ * center + _QUART_ * (below+above)
                  below  = NN(i,j,k)
                  NN(i,j,k) = center
                  center = NN(i,j,k+1)
                  above  = NN(i,j,k+2)
               end do
            end if
#endif
         end if
      end do
   end do
!$OMP END DO

#ifdef SMOOTH_BVF_HORI
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then
            do k=kmax-1,1,-1
               NNc = NN(i,j,k)
               if (az(i+1,j) .ge. 1) then
                  NNe= NN(i+1,j,k)
               else
                  NNe=NNc
               end if
               if (az(i-1,j) .ge. 1) then
                  NNw= NN(i-1,j,k)
               else
                  NNw=NNc
               end if
               if (az(i,j+1) .ge. 1) then
                  NNn= NN(i,j+1,k)
               else
                  NNn=NNc
               end if
               if (az(i,j-1) .ge. 1) then
                  NNs= NN(i,j-1,k)
               else
                  NNs=NNc
               end if
               NN(i,j,k) = (_ONE_/3)*NNc+(_ONE_/6)*(NNe+NNw+NNn+NNs)
            end do
         end if
      end do
   end do
!$OMP END DO
#endif

#ifdef SLICE_MODEL
do k=1,kmax-1
   do i = imin,imax
      if (az(i,2) .ge. 1 ) then
         NN(i,3,k)=NN(i,2,k)
      end if
   end do
end do
#endif

!$OMP END PARALLEL
   call toc(TIM_NN)
#ifdef DEBUG
   write(debug,*) 'Leaving buoyancy_frequency()'
   write(debug,*)
#endif
   return
   end subroutine buoyancy_frequency
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
