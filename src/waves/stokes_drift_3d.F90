#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stokes_drift_3d - calculates layer-integrated Stokes drift
!
! !INTERFACE:
   subroutine stokes_drift_3d(dt,Dveln,hvel,uuEx,vvEx)
!
! !USES:
   use halo_zones     , only: U_TAG,V_TAG
   use domain         , only: imin,imax,jmin,jmax,kmax,au,av
   use pool           , only: flux_center2interface
   use waves          , only: waves_method,WAVES_RS
   use variables_waves, only: waveK,waveE
   use variables_waves, only: SJJ,sinh2kDvelnm1,khab,layerratios
   use variables_waves, only: UStokesCadv,VStokesCadv
   use variables_waves, only: uuStokesC,vvStokesC,uuStokes,vvStokes
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt
   REALTYPE,dimension(I2DFIELD),intent(in)    :: Dveln
   REALTYPE,dimension(I3DFIELD),intent(in)    :: hvel
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: uuEx,vvEx
!
! !DESCRIPTION:
!  Has to be called every 3D timestep after coordinates to take care of
!  adaptively changing layer distribution.
!  stokes_drift_3d is called first for each 3D timestep,
!  so some quantities calculated here are also used later in
!  other routines (sinh2kDcm1,SJJ,khab,layerratios)
!  for RS tendency of stokes transport is added to uuEx here,
!  because then it woll only enter the momentum routines, but not the
!  slow terms
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES
   REALTYPE,dimension(I2DFIELD,0:1) :: sinh2khab
   REALTYPE,dimension(I2DFIELD)     :: hab
   REALTYPE                         :: dtm1
   integer                          :: i,j,k,km,kp
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stokes_drift_3d() # ',Ncall
#endif

   sinh2kDvelnm1 = _ONE_ / sinh(_TWO_*waveK*Dveln)

!  wave-induced kinematic pressure
   SJJ = waveE * waveK * sinh2kDvelnm1

   hab = _ZERO_
   sinh2khab(:,:,0) = _ZERO_
   kp = 0

   do k=1,kmax
      km = kp
      kp = 1 - kp
      hab = hab + hvel(:,:,k)
      khab(:,:,k) = waveK*hab
      sinh2khab  (:,:,kp) = sinh(_TWO_*khab(:,:,k))
      layerratios(:,:,k ) = ( sinh2khab(:,:,kp) - sinh2khab(:,:,km) ) * sinh2kDvelnm1
   end do


!  here begins the actual calculation of the layer-integrated Stokes drift
   dtm1 = _ONE_ / dt


   do k=1,kmax

!     Due to semi-implicit treatment of Eulerian velocity in momentum
!     routines, for RS the tendency of the Stokes transport is added
!     as forcing term.
      if (waves_method .eq. WAVES_RS) then
         do j=jmin,jmax
            do i=imin,imax
               if ( au(i,j).eq.1 .or. au(i,j).eq.2 ) then
                  uuEx(i,j,k) = uuEx(i,j,k) - uuStokes(i,j,k)*dtm1
               end if
               if ( av(i,j).eq.1 .or. av(i,j).eq.2 ) then
                  vvEx(i,j,k) = vvEx(i,j,k) - vvStokes(i,j,k)*dtm1
               end if
            end do
         end do
      end if

      uuStokesC(:,:,k) = layerratios(:,:,k) * UStokesCadv
      call flux_center2interface(U_TAG,uuStokesC(:,:,k),U_TAG,uuStokes(:,:,k))
      call mirror_bdy_2d(uuStokes(:,:,k),U_TAG)

      vvStokesC(:,:,k) = layerratios(:,:,k) * VStokesCadv
      call flux_center2interface(V_TAG,vvStokesC(:,:,k),V_TAG,vvStokes(:,:,k))
      call mirror_bdy_2d(vvStokes(:,:,k),V_TAG)

!     Due to semi-implicit treatment of Eulerian velocity in momentum
!     routines, for RS the tendency of the Stokes transport is added
!     as forcing term.
      if (waves_method .eq. WAVES_RS) then
         do j=jmin,jmax
            do i=imin,imax
               if ( au(i,j).eq.1 .or. au(i,j).eq.2 ) then
                  uuEx(i,j,k) = uuEx(i,j,k) + uuStokes(i,j,k)*dtm1
               end if
               if ( av(i,j).eq.1 .or. av(i,j).eq.2 ) then
                  vvEx(i,j,k) = vvEx(i,j,k) + vvStokes(i,j,k)*dtm1
               end if
            end do
         end do
      end if

   end do


#ifdef DEBUG
   write(debug,*) 'Leaving stokes_drift_3d()'
   write(debug,*)
#endif
   return
   end subroutine stokes_drift_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
