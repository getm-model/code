#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pool
!
! !INTERFACE:
   module pool
!
! !DESCRIPTION:
!
!
! !USES:
   use domain
   use exceptions

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   interface

      subroutine deformation_rates(U,V,DU,DV,kwe,dudxC,dudxV,dudxU,         &
                                                 dvdyC,dvdyU,dvdyV,         &
                                                 dudyX,dvdxX,shearX,        &
                                                 dvdxU,shearU,dudyV,shearV)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
         logical,intent(in),optional                       :: kwe !keyword-enforcer
         REALTYPE,dimension(E2DFIELD),target,intent(out),optional :: dudxC,dvdyC
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxV,dudxU
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dvdyU,dvdyV
         REALTYPE,dimension(:,:),pointer,intent(out),optional :: dudyX,dvdxX
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: shearX
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dvdxU,shearU
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyV,shearV
      end subroutine deformation_rates

      subroutine flux_center2interface(tagc,fluxc,ttag,fluxi)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         integer,intent(in)                       :: tagc,ttag
         REALTYPE,dimension(E2DFIELD),intent(in)  :: fluxc
         REALTYPE,dimension(E2DFIELD),intent(out) :: fluxi
      end subroutine flux_center2interface

   end interface
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   end module pool

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
