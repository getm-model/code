#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: radiation_stress_3d - calculates layer-integrated Radiation Stress
!
! !INTERFACE:
   subroutine radiation_stress_3d(Dveln,hvel,uuEx,vvEx)
!
! !USES:
   use domain         , only: imin,imax,jmin,jmax,kmax,az
   use waves          , only: kD_max
   use variables_waves, only: waveK,waveE
   use variables_waves, only: SJJ
   use variables_waves, only: kDveln,mask_kDveln
   use variables_waves, only: khab,layerratios
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(I2DFIELD),intent(in)    :: Dveln
   REALTYPE,dimension(I3DFIELD),intent(in)    :: hvel
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: uuEx,vvEx
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES
   REALTYPE,dimension(I2DFIELD,0:1) :: sinhkhab2
   REALTYPE,dimension(I2DFIELD)     :: sinhkDvelnm2,layerratios2,SE,gradterms
   integer                          :: i,j,k,km,kp
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'radiation_stress_3d() # ',Ncall
#endif

   sinhkDvelnm2 = (_ONE_ / sinh(kDveln))**2

   sinhkhab2(:,:,0) = _ZERO_
   kp = 0


   do k=1,kmax

      km = kp
      kp = 1 - kp

      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if ( az(i,j) .gt. 0 ) then
               if ( mask_kDveln(i,j) ) then
                  sinhkhab2   (i,j,kp) = sinh(khab(i,j,k))**2
                  layerratios2(i,j) = ( sinhkhab2(i,j,kp) - sinhkhab2(i,j,km) ) * sinhkDvelnm2(i,j)
               else
                  layerratios2(i,j) = layerratios(i,j,k)
               end if
            end if
         end do
      end do

      SE  = _HALF_*waveE*layerratios(:,:,k) + hvel(:,:,k)*SJJ
      gradterms = SE - _HALF_*waveE*layerratios2

      call rs_force(SE,gradterms,uuEx(:,:,k),vvEx(:,:,k))

   end do


#ifdef DEBUG
   write(debug,*) 'Leaving radiation_stress_3d()'
   write(debug,*)
#endif
   return
   end subroutine radiation_stress_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
