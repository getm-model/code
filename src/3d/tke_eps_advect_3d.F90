#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tke_eps_advect_3d - 3D turbulence advection
!
! !INTERFACE:
   subroutine tke_eps_advect_3d()
!
! !DESCRIPTION:
!
! This routine carries out advection of the prognostic turbulence quantities
! {\tt tke} (turbuent kinetic energy, $k$) and {\tt eps} (lenght scale related
! turbulence quantity, e.g.\ dissipation rate of $k$, $\varepsilon$, or
! turbulent frequency, $\omega=\varepsilon/k$. Here, the TVD advection
! schemes are used which are also used for the momentum advection.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,az
   use m3d, only: vel_adv_split,vel_hor_adv,vel_ver_adv
   use variables_3d, only: tke,eps,dt,uu,vv,ww,hun,hvn,ho,hn
   use variables_3d, only: uuadv,vvadv,wwadv,hoadv,hnadv,huadv,hvadv
   use advection_3d, only: do_advection_3d
   use halo_zones, only: H_TAG
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tke_eps_advect_3d() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

   do k=1,kmax-1
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            uuadv(i,j,k) = _HALF_*( uu(i,j,k) + uu(i,j,k+1) )
            vvadv(i,j,k) = _HALF_*( vv(i,j,k) + vv(i,j,k+1) )
            wwadv(i,j,k) = _HALF_*( ww(i,j,k) + ww(i,j,k+1) )
            hoadv(i,j,k) = _HALF_*( ho(i,j,k) + ho(i,j,k+1) )
            hnadv(i,j,k) = _HALF_*( hn(i,j,k) + hn(i,j,k+1) )
!           Note (KK): hun only valid until imax+1
!                      therefore huadv only valid until imax+1
            huadv(i,j,k) = _HALF_*( hun(i,j,k) + hun(i,j,k+1) )
!           Note (KK): hvn only valid until jmax+1
!                      therefore hvadv only valid until jmax+1
            hvadv(i,j,k) = _HALF_*( hvn(i,j,k) + hvn(i,j,k+1) )
         end do
      end do
!$OMP END DO NOWAIT
   end do

!  only needed to allow for advection of k=kmax
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
!        KK-TODO: does _HALF_ transports make sense?
         uuadv(i,j,kmax) = _HALF_*uu(i,j,kmax)
         vvadv(i,j,kmax) = _HALF_*vv(i,j,kmax)
         wwadv(i,j,kmax) = _HALF_*ww(i,j,kmax)
         hoadv(i,j,kmax) = _HALF_*ho(i,j,kmax)
         hnadv(i,j,kmax) = _HALF_*hn(i,j,kmax)
         huadv(i,j,kmax) = _HALF_*hun(i,j,kmax)
         hvadv(i,j,kmax) = _HALF_*hvn(i,j,kmax)
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

   call do_advection_3d(dt,tke,uuadv,vvadv,wwadv,huadv,hvadv,hoadv,hnadv,   &
                        vel_hor_adv,vel_ver_adv,vel_adv_split,_ZERO_,H_TAG)

   call do_advection_3d(dt,eps,uuadv,vvadv,wwadv,huadv,hvadv,hoadv,hnadv,   &
                        vel_hor_adv,vel_ver_adv,vel_adv_split,_ZERO_,H_TAG)

#ifdef DEBUG
   write(debug,*) 'Leaving tke_eps_advect_3d()'
   write(debug,*)
#endif
   return
   end subroutine tke_eps_advect_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
