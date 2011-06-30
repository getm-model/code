!$Id: calc_uvex.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: calc_uvex
!
! !INTERFACE:
   subroutine calc_uvex(U,V,D,DU,DV)
!
! !DESCRIPTION:
!
! !USES:
   use m2d, only: Am_method,An_method
   use variables_2d, only: UEx,VEx
   use variables_2d, only: dudxC,dvdyC,shearX
   use variables_2d, only: dudxU,dvdyV,shearU
#ifdef _LES_
   use les, only: do_les_2d
   use variables_les, only: Am_2d,AmX_2d
#endif
!  needed by interface headers!
   use domain, only: imin,imax,jmin,jmax

   IMPLICIT NONE

#include "deformation_rates.h"
#include "uv_diffusion.h"

!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: U,V,D,DU,DV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'calc_uvex() # ',Ncall
#endif

      call uv_advect(U,V,DU,DV)
      select case(Am_method)
         case(0)
            if (An_method .gt. 0) then
               call uv_diffusion(0,An_method,UEx,VEx,U=U,V=V)
            end if
         case(1)
            call deformation_rates(U,V,DU,DV,dudxC,dvdyC,shearX=shearX)
            call uv_diffusion(1,An_method,UEx,VEx,U=U,V=V,           &
                              dudxC=dudxC,dvdyC=dvdyC,shearX=shearX, &
                              D=D,DU=DU,DV=DV)
#ifdef _LES_
         case(2)
            call deformation_rates(U,V,DU,DV,dudxC,dvdyC,shearX=shearX, &
                                   dudxU=dudxU,dvdyV=dvdyV,shearU=shearU)
            call do_les_2d(dudxC,dudxU,dvdyC,dvdyV,shearX,shearU)
            call uv_diffusion(2,An_method,UEx,VEx,U=U,V=V,           &
                              dudxC=dudxC,dvdyC=dvdyC,shearX=shearX, &
                              D=D,DU=DU,DV=DV,Am=Am_2d,AmX=AmX_2d)
#endif
      end select

#ifdef DEBUG
   write(debug,*) 'Leaving calc_uvex()'
   write(debug,*)
#endif
   return
   end subroutine calc_uvex

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
