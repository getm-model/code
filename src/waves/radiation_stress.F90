#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: radiation_stress - calculates depth-integrated Radiation Stress
!
! !INTERFACE:
   subroutine radiation_stress(Dvel,UEx,VEx)
!
! !USES:
   use domain         , only: imin,imax,jmin,jmax
   use variables_waves, only: waveE,SJ
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)    :: Dvel
!
! !INPUT/OUTPUT PARAMETERS:
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
   REALTYPE,dimension(E2DFIELD) :: DSJ,SE
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'radiation_stress() # ',Ncall
#endif

   DSJ = Dvel * SJ
   SE  = _HALF_*waveE + DSJ

   call rs_force(SE,DSJ,UEx,VEx)

#ifdef DEBUG
   write(debug,*) 'Leaving radiation_stress()'
   write(debug,*)
#endif
   return
   end subroutine radiation_stress
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
