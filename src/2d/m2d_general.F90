#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m2d_general - interfaces for 2d subroutines
!
! !INTERFACE:
   module m2d_general
!
! !DESCRIPTION:
!
! !USE:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine uv_diffusion(An_method,UEx,VEx,U,V,D,DU,DV,hsd_u,hsd_v)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         integer,intent(in)                                :: An_method
         REALTYPE,dimension(E2DFIELD),intent(in),optional  :: U,V,D,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(inout)        :: UEx,VEx
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: hsd_u,hsd_v
      end subroutine uv_diffusion

      subroutine calc_uvex(An_method,U,V,D,DU,DV)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         integer,intent(in)                      :: An_method
         REALTYPE,dimension(E2DFIELD),intent(in) :: U,V,D,DU,DV
      end subroutine calc_uvex

      subroutine bottom_friction(U,V,DU,DV,ru,rv,zub,zvb)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(out)          :: ru,rv
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: zub,zvb
      end subroutine bottom_friction
   end interface

!-----------------------------------------------------------------------

   end module m2d_general

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
