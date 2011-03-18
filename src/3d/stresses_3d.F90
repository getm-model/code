#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stresses_3d - bottom and surface stresses\label{sec-stresses-3d}
!
! !INTERFACE:
   subroutine stresses_3d
!
! !DESCRIPTION:
!
! As preparation of the call to {\tt do\_turbulence} in the routine {\tt gotm},
! see section \ref{sec-gotm}, the normalised surface and bottom stresses,
! $\tau_s/\rho_0$ (variable {\tt taus}) and $\tau_b/\rho_0$ (variable {\tt
! taub}), respectively, are calculated and interpolated to the T-points.
! Input parameters to this routine are {\tt rru} and {tt rrv}, which
! contain $r\sqrt{u^2+v^2}$ for the U- and V-points, respectively.
! The modules of the surface and bottom stress vectors are calculated
! then by means of taking the square root of the sum of the squares of
! the stess components. In a similar way also the $x$- and $y$-components
! of the bottom stress are computed for output.
!
! !USES:
   use parameters, only: rho_0
   use domain, only: az,au,av,imin,imax,jmin,jmax
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,rru,rrv
   use variables_3d, only: taus,taubx,tauby,taub
   use meteo, only: tausx,tausy
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
   use getm_timers, only: tic, toc, TIM_STRESSES3D, TIM_STRESSES3DH
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,ku1,ku2,kv1,kv2
   REALTYPE                  :: rho_0i
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stresses_3d() # ',Ncall
#endif
   call tic(TIM_STRESSES3D)

   rho_0i=_ONE_/rho_0

!  we need to know rru and rrv in the halos as well
   call tic(TIM_STRESSES3DH)
   call update_2d_halo(rru,rru,au,imin,jmin,imax,jmax,10)
   call wait_halo(10)
   call update_2d_halo(rrv,rrv,av,imin,jmin,imax,jmax,10)
   call wait_halo(10)
   call toc(TIM_STRESSES3DH)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ku1,ku2,kv1,kv2)

!  x-component of bottom momentum flux at U-points
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax

         k          = kumin(i,j)                      ! bottom index
         taubx(i,j) = - uu(i,j,k)/hun(i,j,k)*rru(i,j) ! momentum flux

      enddo
   enddo
!$OMP END DO NOWAIT

!  y-component of bottom momentum flux at V-points
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax

         k          = kvmin(i,j)                      ! bottom index
         tauby(i,j) = - vv(i,j,k)/hvn(i,j,k)*rrv(i,j) ! momentum flux

      enddo
   enddo
!$OMP END DO

!  stress magnitude
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax

!        lower indices at U- and V-points
         ku1=kumin(i-1,j  )
         ku2=kumin(i  ,j  )
         kv1=kvmin(i  ,j-1)
         kv2=kvmin(i  ,j  )


!        total bottom stress at T-points
         taub(i,j)=sqrt(_HALF_*(                                      &
                   (uu(i-1,j  ,ku1)/hun(i-1,j  ,ku1)*rru(i-1,j  ))**2 &
                  +(uu(i  ,j  ,ku2)/hun(i  ,j  ,ku2)*rru(i  ,j  ))**2 &
                  +(vv(i  ,j-1,kv1)/hvn(i  ,j-1,kv1)*rrv(i  ,j-1))**2 &
                  +(vv(i  ,j  ,kv2)/hvn(i  ,j  ,kv2)*rrv(i  ,j  ))**2))

!        total surface stress at T-points
         taus(i,j)=rho_0i*sqrt(_HALF_*(                 &
                        tausx(i,j)**2+tausx(i-1,j)**2   &
                      + tausy(i,j)**2+tausy(i,j-1)**2) )

      end do
   end do
!$OMP END DO
!$OMP END PARALLEL

   call toc(TIM_STRESSES3D)
#ifdef DEBUG
   write(debug,*) 'Leaving stresses_3d()'
   write(debug,*)
#endif
   return
   end subroutine stresses_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
