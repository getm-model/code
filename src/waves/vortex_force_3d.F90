#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: vortex_force_3d - layer-integrated VF
!
! !INTERFACE:
   subroutine vortex_force_3d(uuEuler,vvEuler,hun,hvn,uuEx,vvEx)
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones     , only: U_TAG,V_TAG
   use domain         , only: imin,imax,jmin,jmax,kmax,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain         , only: dxu,dyv
#else
   use domain         , only: dx,dy
#endif
   use pool           , only: deformation_rates,flux_center2interface
   use variables_waves, only: SJJ
   use variables_waves, only: uuStokesC,vvStokesC,uuStokes,vvStokes
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(in)    :: uuEuler,vvEuler,hun,hvn
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: uuEx,vvEx
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(I2DFIELD) :: dJdx,dJdy,dudxU,dvdyV,dvdxU,dudyV
   REALTYPE,dimension(I2DFIELD) :: work2d
   integer                      :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Wcall = 0
   Wcall = Wcall+1
   write(debug,*) 'vortex_force_3d() # ',Wcall
#endif

!  wave-induced pressure gradient at U-points
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            dJdx(i,j) = ( SJJ(i+1,j) - SJJ(i,j) ) / DXU
         end if
      end do
   end do

!  wave-induced pressure gradient at V-points
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            dJdy(i,j) = ( SJJ(i,j+1) - SJJ(i,j) ) / DYV
         end if
      end do
   end do


   do k=1,kmax

!     KK-TODO: use of already calculated deformation rates...
      call deformation_rates(uuEuler(:,:,k),vvEuler(:,:,k),hun(:,:,k),hvn(:,:,k), &
                             dudxU=dudxU,dvdyV=dvdyV,dvdxU=dvdxU,dudyV=dudyV)
!     KK-TODO: correct from layer-orientated gradients to z-coordinates

!     layer-integrated Stokes drift in y-direction at U-point
      call flux_center2interface(V_TAG,vvStokesC(:,:,k),U_TAG,work2d)

      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j).eq.1 .or. au(i,j).eq.2) then
               uuEx(i,j,k) =   uuEx(i,j,k)                &
                             - uuStokes(i,j,k)*dudxU(i,j) &
                             - work2d  (i,j  )*dvdxU(i,j) &
                             + hun(i,j,k)*dJdx(i,j)
                                
            end if
         end do
      end do

!     layer-integrated Stokes drift in x-direction at V-point
      call flux_center2interface(U_TAG,uuStokesC(:,:,k),V_TAG,work2d)

      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j).eq.1 .or. av(i,j).eq.2) then
               vvEx(i,j,k) =   vvEx(i,j,k)                &
                             - vvStokes(i,j,k)*dvdyV(i,j) &
                             - work2d  (i,j  )*dudyV(i,j) &
                             + hvn(i,j,k)*dJdy(i,j)
            end if
         end do
      end do

!     KK-TODO: add dissipation terms at surface and bottom

   end do

   end subroutine vortex_force_3d
!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2013 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
