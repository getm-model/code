#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: vortex_force - depth-integrated VF
!
! !INTERFACE:
   subroutine vortex_force(UEuler,VEuler,DU,DV,UEx,VEx)
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones     , only: U_TAG,V_TAG
   use domain         , only: imin,imax,jmin,jmax,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain         , only: dxu,dyv
#else
   use domain         , only: dx,dy
#endif
   use pool           , only: deformation_rates,flux_center2interface
   use variables_waves, only: SJ
   use variables_waves, only: UStokesC,VStokesC,UStokes,VStokes
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)    :: UEuler,VEuler,DU,DV
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: UEx,VEx
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD) :: dJdx,dJdy,dudxU,dvdyV,dvdxU,dudyV
   REALTYPE,dimension(E2DFIELD) :: work2d
   integer                      :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Wcall = 0
   Wcall = Wcall+1
   write(debug,*) 'vortex_force() # ',Wcall
#endif

!  wave-induced pressure gradient at U-points
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            dJdx(i,j) = ( SJ(i+1,j) - SJ(i,j) ) / DXU
         end if
      end do
   end do

!  wave-induced pressure gradient at V-points
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            dJdy(i,j) = ( SJ(i,j+1) - SJ(i,j) ) / DYV
         end if
      end do
   end do

!  KK-TODO: use of already calculated deformation rates...
   call deformation_rates(UEuler,VEuler,DU,DV,                             &
                          dudxU=dudxU,dvdyV=dvdyV,dvdxU=dvdxU,dudyV=dudyV)

!  depth-integrated Stokes drift in y-direction at U-point
   call flux_center2interface(V_TAG,VStokesC,U_TAG,work2d)

   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            UEx(i,j) =   UEx(i,j)                &
                       - UStokes(i,j)*dudxU(i,j) &
                       - work2d (i,j)*dvdxU(i,j) &
                       + DU(i,j)*dJdx(i,j)
                             
         end if
      end do
   end do

!  depth-integrated Stokes drift in x-direction at V-point
   call flux_center2interface(U_TAG,UStokesC,V_TAG,work2d)

   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            VEx(i,j) =   VEx(i,j)                &
                       - VStokes(i,j)*dvdyV(i,j) &
                       - work2d (i,j)*dudyV(i,j) &
                       + DV(i,j)*dJdy(i,j)
         end if
      end do
   end do

!  KK-TODO: add dissipation terms at surface and bottom

   end subroutine vortex_force
!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2013 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
