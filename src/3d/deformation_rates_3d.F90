!$Id: deformation_rates_3d.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: deformation_rates_3d - \label{deformation_rates_3d}
!
! !INTERFACE:
   subroutine deformation_rates_3d()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: kmax
   use m3d, only: deformUV
   use variables_3d, only: uu,vv,hun,hvn
   use variables_3d, only: dudxC_3d,dvdyC_3d,shearX_3d
   use variables_3d, only: dudxU_3d,dvdyV_3d,shearU_3d

   IMPLICIT NONE

#include "../2d/deformation_rates.h"
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer                                           :: k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'deformation_rates_3d() # ',Ncall
#endif

!  KK-TODO: timer
   do k=1,kmax
      if (deformUV) then
         call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k),     &
                                dudxC_3d(:,:,k),dvdyC_3d(:,:,k),               &
                                dudxU=dudxU_3d(:,:,k),dvdyV=dvdyV_3d(:,:,k),   &
                                shearX=shearX_3d(:,:,k),shearU=shearU_3d(:,:,k))
      else
         call deformation_rates(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k), &
                                dudxC_3d(:,:,k),dvdyC_3d(:,:,k),           &
                                shearX=shearX_3d(:,:,k))
      end if
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving deformation_rates_3d()'
   write(debug,*)
#endif
   return
   end subroutine deformation_rates_3d

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
