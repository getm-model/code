#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  REALTYPE function adv_tvd_limiter -
!
! !INTERFACE:
   REALTYPE function adv_tvd_limiter(scheme,cfl,slope)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! !USES:
   use advection, only: P2,SUPERBEE,MUSCL,P2_PDM
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)  :: scheme
   REALTYPE,intent(in) :: cfl,slope
!
! !LOCAL VARIABLES:
   REALTYPE           :: x,Phi
   REALTYPE,parameter :: one6th=_ONE_/6
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (scheme)
      case ((P2),(P2_PDM))
         x = one6th*(_ONE_-_TWO_*cfl)
         Phi = (_HALF_+x) + (_HALF_-x)*slope
         if (scheme.eq.P2) then
            adv_tvd_limiter = Phi
         else
            adv_tvd_limiter = max(_ZERO_,min(Phi,_TWO_/(_ONE_-cfl),_TWO_*slope/(cfl+1.d-10)))
         end if
      case (SUPERBEE)
         adv_tvd_limiter = max(_ZERO_,min(_ONE_,_TWO_*slope),min(slope,_TWO_))
      case (MUSCL)
         adv_tvd_limiter = max(_ZERO_,min(_TWO_,_TWO_*slope,_HALF_*(_ONE_+slope)))
      case default
         stop 'adv_tvd_limiter: invalid scheme'
   end select

   return
   end function adv_tvd_limiter
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
