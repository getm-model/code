#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stokes_drift - calculates depth-integrated Stokes drift
!
! !INTERFACE:
   subroutine stokes_drift(dtm,Dvel,UEx,VEx)
!
! !USES:
   use halo_zones     , only: U_TAG,V_TAG
   use domain         , only: imin,imax,jmin,jmax,au,av
   use pool           , only: flux_center2interface
   use waves          , only: new_StokesC,waves_method,WAVES_RS
   use variables_waves, only: waveE,waveK,SJ
   use variables_waves, only: UStokesC,VStokesC,UStokes,VStokes
   IMPLICIT NONE
!
! !INPUT VARIABLES:
   REALTYPE,intent(in)                        :: dtm
   REALTYPE,dimension(E2DFIELD),intent(in)    :: Dvel
!
! !INPUT/OUTPUT VARIABLES:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: UEx,VEx
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES
   REALTYPE :: dtmm1
   integer  :: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stokes_drift() # ',Ncall
#endif

!  stokes_drift is called first for each timestep,
!  so some quantities calculated here are also used later in
!  other routines (SJ)

!  wave-induced kinematic pressure
   SJ = waveE * waveK / sinh(_TWO_*waveK*Dvel)


!  here begins the actual calculation of the depth-integrated Stokes drift
   if (new_StokesC) then

      new_StokesC = .false.

!     Due to semi-implicit treatment of Eulerian velocity in momentum
!     routines, for RS the tendency of the Stokes transport is added
!     as forcing term.
      if (waves_method .eq. WAVES_RS) then
         dtmm1 = _ONE_ / dtm
         do j=jmin,jmax
            do i=imin,imax
               if ( au(i,j).eq.1 .or. au(i,j).eq.2 ) then
                  UEx(i,j) = UEx(i,j) - UStokes(i,j)*dtmm1
               end if
               if ( av(i,j).eq.1 .or. av(i,j).eq.2 ) then
                  VEx(i,j) = VEx(i,j) - VStokes(i,j)*dtmm1
               end if
            end do
         end do
      end if

!     average to U-point
      call flux_center2interface(U_TAG,UStokesC,U_TAG,UStokes)
      call mirror_bdy_2d(UStokes,U_TAG)

!     average to V-point
      call flux_center2interface(V_TAG,VStokesC,V_TAG,VStokes)
      call mirror_bdy_2d(VStokes,V_TAG)

!     Due to semi-implicit treatment of Eulerian velocity in momentum
!     routines, for RS the tendency of the Stokes transport is added
!     as forcing term.
      if (waves_method .eq. WAVES_RS) then
         do j=jmin,jmax
            do i=imin,imax
               if ( au(i,j).eq.1 .or. au(i,j).eq.2 ) then
                  UEx(i,j) = UEx(i,j) + UStokes(i,j)*dtmm1
               end if
               if ( av(i,j).eq.1 .or. av(i,j).eq.2 ) then
                  VEx(i,j) = VEx(i,j) + VStokes(i,j)*dtmm1
               end if
            end do
         end do
      end if

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving stokes_drift()'
   write(debug,*)
#endif
   return
   end subroutine stokes_drift
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
