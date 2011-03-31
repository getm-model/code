#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tke_eps_advect_3d - 3D turbulence advection
!
! !INTERFACE:
   subroutine tke_eps_advect_3d(hor_adv,ver_adv,adv_split)
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
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,kumin,kvmin,uu,vv,ww,hun,hvn,ho,hn,uuEx,vvEx
#ifdef UV_TVD
   use variables_3d, only: uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv
#endif
   use variables_3d, only: tke,eps
   use advection_3d, only: do_advection_3d
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG,V_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: hor_adv,ver_adv,adv_split
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
#ifdef UV_TVD
   REALTYPE                  :: dxuadv(I2DFIELD)
   REALTYPE                  :: dxvadv(I2DFIELD)
   REALTYPE                  :: dyuadv(I2DFIELD)
   REALTYPE                  :: dyvadv(I2DFIELD)
   REALTYPE                  :: area_inv(I2DFIELD)
   REALTYPE                  :: AH=_ZERO_
   REALTYPE                  :: dxdyi
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tke_eps_advect_3d() # ',Ncall
#endif

#ifndef UV_TVD
   STDERR 'Do not use tke_eps_advect_3d() without the compiler option UV_TVD'
   stop 'tke_eps_advect_3d()'
#else

#if defined(SPHERICAL) || defined(CURVILINEAR)
#else
   dxdyi = _ONE_/(dx*dy)
#endif

! Here begins dimensional split advection for tke and eps
   do k=1,kmax-1
     do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            uadv(i,j,k)=_HALF_*(uu(i,j,k)+uu(i,j,k+1))
            vadv(i,j,k)=_HALF_*(vv(i,j,k)+vv(i,j,k+1))
            wadv(i,j,k)=_HALF_*(ww(i,j,k)+ww(i,j,k+1))
            huadv(i,j,k)=_HALF_*(hun(i,j,k)+hun(i,j,k+1))
            hvadv(i,j,k)=_HALF_*(hvn(i,j,k)+hvn(i,j,k+1))
            hoadv(i,j,k)=_HALF_*(ho(i,j,k)+ho(i,j,k+1))
            hnadv(i,j,k)=_HALF_*(hn(i,j,k)+hn(i,j,k+1))
         end do
      end do
   end do
   k=kmax  ! Only needed to allow for advection of k=kmax
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         uadv(i,j,k)=_HALF_*(uu(i,j,k))
         vadv(i,j,k)=_HALF_*(vv(i,j,k))
         wadv(i,j,k)=_HALF_*(ww(i,j,k))
         huadv(i,j,k)=_HALF_*(hun(i,j,k))
         hvadv(i,j,k)=_HALF_*(hvn(i,j,k))
         hoadv(i,j,k)=_HALF_*(ho(i,j,k))
         hnadv(i,j,k)=_HALF_*(hn(i,j,k))
      end do
   end do

   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
#if defined(SPHERICAL) || defined(CURVILINEAR)
         dxuadv(i,j)=dxu(i,j)
         dxvadv(i,j)=dxv(i,j)
         dyuadv(i,j)=dyu(i,j)
         dyvadv(i,j)=dyv(i,j)
         area_inv(i,j)=arcd1(i,j)
#else
         dxuadv(i,j)=dx
         dxvadv(i,j)=dx
         dyuadv(i,j)=dy
         dyvadv(i,j)=dy
         area_inv(i,j)=dxdyi
#endif
      end do
   end do

   call do_advection_3d(dt,tke,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv, &
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,          &
                        az,au,av,hor_adv,ver_adv,adv_split,AH)

   call do_advection_3d(dt,eps,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv, &
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,          &
                        az,au,av,hor_adv,ver_adv,adv_split,AH)

#endif

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
