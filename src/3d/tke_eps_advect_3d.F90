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
   use domain, only: imin,imax,jmin,jmax,kmax,az,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxv,dyu
#else
   use domain, only: dx,dy
#endif
   use m3d, only: vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver
   use variables_3d, only: tke,eps,dt,uu,vv,ww,hun,hvn,ho,hn
   use advection, only: J7
   use advection_3d, only: do_advection_3d,W_TAG
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG
   use turbulence, only: k_min,eps_min
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                      :: i,j,k
   REALTYPE,dimension(I3DFIELD) :: uuadv,vvadv,wwadv,hoadv,hnadv,huadv,hvadv
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tke_eps_advect_3d() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif

!  Note (KK): in the present implementation it is guaranteed, that 
!             the turbulent quantities at k=0|kmax are not changed
!             due to advection
!   uuadv(:,:,0) = _ZERO_ ! not used in the present implementation
!   vvadv(:,:,0) = _ZERO_ ! not used in the present implementation
   wwadv(:,:,0) = _HALF_*ww(:,:,1)
   hoadv(:,:,0) = _HALF_*ho(:,:,1)
!   hnadv(:,:,0) = _HALF_*hn(:,:,1) ! not used in the present implementation
!   huadv(:,:,0) = _HALF_*hun(:,:,1) ! not used in the present implementation
!   hvadv(:,:,0) = _HALF_*hvn(:,:,1) ! not used in the present implementation

   uuadv(:,:,kmax) = _ZERO_
   vvadv(:,:,kmax) = _ZERO_
!   wwadv(:,:,kmax) = _ZERO_ ! not used in the present implementation
   hoadv(:,:,kmax) = _HALF_*ho(:,:,kmax)
   hnadv(:,:,kmax) = _HALF_*hn(:,:,kmax)
   huadv(:,:,kmax) = _HALF_*hun(:,:,kmax)
   hvadv(:,:,kmax) = _HALF_*hvn(:,:,kmax)


!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k)

   do k=1,kmax-1
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
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
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT
   end do

!$OMP END PARALLEL

   if (vel3d_adv_hor .eq. J7) then
#ifdef SLICE_MODEL
      uuadv(:,j+1,:) = uuadv(:,j,:)
#endif
      do k=1,kmax-1
!        OMP-NOTE: these loops cannot be threaded!
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO-1
#endif
            do i=imin-HALO,imax+HALO-1
               if (ax(i,j) .eq. 0) then
                  uuadv(i,j,k) = _ZERO_
                  vvadv(i,j,k) = _ZERO_
               else
!                 Note (KK): due to mirroring open bdys are already included here
                  uuadv(i,j,k) = _HALF_*( uuadv(i,j,k)*DYU + uuadv(i,j+1,k)*DYUJP1 )
                  vvadv(i,j,k) = _HALF_*( vvadv(i,j,k)*DXV + vvadv(i+1,j,k)*DXVIP1 )
               end if
!              KK-TODO: do we have to average h[o|n]adv as well?
            end do
#ifndef SLICE_MODEL
         end do
#endif
      end do
   end if

!  Note (KK): tke and eps are not halo-updated yet
   call update_3d_halo(tke,tke,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   call update_3d_halo(eps,eps,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   call do_advection_3d(dt,tke,uuadv,vvadv,wwadv,huadv,hvadv,hoadv,hnadv,         &
                        vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,W_TAG)

   call do_advection_3d(dt,eps,uuadv,vvadv,wwadv,huadv,hvadv,hoadv,hnadv,         &
                        vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,W_TAG)

   tke = max(k_min,tke)
   eps = max(eps_min,eps)

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
