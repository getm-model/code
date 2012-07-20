#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: shear_frequency - calculates the shear frequency
! \label{sec-shear-frequency}
!
! !INTERFACE:
   subroutine shear_frequency()
!
! !DESCRIPTION:
!
! Here, the shear frequency squared,
! $M^2=\left(\partial_z u\right)^2+\left(\partial_z v\right)^2$,
! is calculated.
! Two alternative methods are coded.
! The straight-forward method which is explained first,
! has the disadvantage of generating numerical instabilities.
! The straight-forward way for calculating $M^2$ is as follows:
!
! \begin{equation}\label{ShearSquaredOld}
! \begin{array}{l}
! \displaystyle
! (M^2)_{i,j,k}\approx
! \frac12
! \Bigg(\left(\frac{u_{i,j,k+1}-u_{i,j,k}}
! {\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}\right)^2
! +
! \left(\frac{u_{i-1,j,k+1}-u_{i-1,j,k}}
! {\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}\right)^2
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \left(\frac{v_{i,j,k+1}-v_{i,j,k}}
! {\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}\right)^2
! +
! \left(\frac{v_{i,j-1,k+1}-v_{i,j-1,k}}
! {\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}\right)^2
! \Bigg)
! \end{array}
! \end{equation}
!
! \cite{BURCHARD01c} developed a new scheme, which guarantees that
! the mean kinetic energy which is dissipated from the mean flow
! equals the shear production of turbulent kinetic energy. Therefore,
! this scheme should be numerically more stable than (\ref{ShearSquaredOld}):
!
! \begin{equation}\label{ShearSquaredNew}
! \begin{array}{l}
! \displaystyle
! (M^2)_{i,j,k}\approx
! \frac12
! \Bigg(\frac{\frac12(\nu_{i,j,k}+\nu_{i+1,j,k})
! (u_{i,j,k+1}-u_{i,j,k})^2}{\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i-1,j,k}+\nu_{i,j,k})
! (u_{i-1,j,k+1}-u_{i-1,j,k})^2}{\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i,j,k}+\nu_{i,j+1,k})
! (v_{i,j,k+1}-v_{i,j,k})^2}{\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i,j-1,k}+\nu_{i,j,k})
! (v_{i,j-1,k+1}-v_{i,j-1,k})^2}{\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}
! \Bigg)
! \\ \\ \qquad\qquad\quad
! \displaystyle
! \cdot
! \left(\frac12\left(h^c_{i,j,k}+h^c_{i,j,k+1}\right)\nu_{i,j,k}\right)^{-1}
! \end{array}
! \end{equation}
!
! In some cases, together with the straight-forward discretisation
! of the buoyancy frequency squared, (\ref{Nstraight}), this
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
! These stability issues need to be further investigated in the future.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,au,av,az
   use variables_3d, only: hn,uu,hun,vv,hvn,SS,num
   use getm_timers, only: tic, toc, TIM_SS
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'shear_frequency() # ',Ncall
#endif
   call tic(TIM_SS)

! BJB-TODO: There are a lot of single- (or lower-) precision constants
!   in this routine. Enough that it significantly influences the result
!   if they are correctly converted to full (double) precision.
!   A conversion should be done, but I will not do it as part of
!   the OMP-implementation as it muddles the comparisons of old and
!   new results. BJB 2009-09-22

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
!  Prandtl frequency
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1 ) then
            SS(i,j,kmax) = _ZERO_
            do k=1,kmax-1
! This is an older version which we should keep here.
#ifndef NEW_SS
              SS(i,j,k)=_HALF_* (                                             &
                   ( (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))          &
                   /(_HALF_*(hun(i,j,k+1)+hun(i,j,k))) )**2                   &
                +  ( (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k))  &
                   /(_HALF_*(hun(i-1,j,k+1)+hun(i-1,j,k))) )**2               &
                +  ( (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))          &
                   /(_HALF_*(hvn(i,j,k+1)+hvn(i,j,k))) )**2                   &
                +  ( (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k))  &
                   /(_HALF_*(hvn(i,j-1,k+1)+hvn(i,j-1,k))) )**2               &
                            )
#else
! This version should better conserve energy.
              SS(i,j,k)=_HALF_* (                                             &
                   (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))**2         &
                   /(_HALF_*(hun(i,j,k+1)+hun(i,j,k)))                        &
                    *_HALF_*(num(i,j,k)+num(i+1,j,k))                         &
               +  (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k))**2  &
                  /(_HALF_*(hun(i-1,j,k+1)+hun(i-1,j,k)))                     &
                   *_HALF_*(num(i-1,j,k)+num(i,j,k))                          &
                +  (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))**2         &
                   /(_HALF_*(hvn(i,j,k+1)+hvn(i,j,k)))                        &
                    *_HALF_*(num(i,j,k)+num(i,j+1,k))                         &
                +  (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k))**2 &
                   /(_HALF_*(hvn(i,j-1,k+1)+hvn(i,j-1,k)))                    &
                    *_HALF_*(num(i,j-1,k)+num(i,j,k))                         &
                            )/(_HALF_*(hn(i,j,k)+hn(i,j,k+1)))/num(i,j,k)
#endif
            end do
         end if
      end do
   end do
!$OMP END DO

#ifdef SLICE_MODEL
do k=1,kmax-1
   do i = imin,imax
      if (az(i,2) .ge. 1 ) then
         SS(i,3,k)=SS(i,2,k)
      end if
   end do
end do
#endif


!$OMP END PARALLEL
   call toc(TIM_SS)
#ifdef DEBUG
   write(debug,*) 'Leaving shear_frequency()'
   write(debug,*)
#endif
   return
   end subroutine shear_frequency
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
