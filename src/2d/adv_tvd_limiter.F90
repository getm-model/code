#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  REALTYPE function adv_tvd_limiter -
!
! !INTERFACE:
   REALTYPE function adv_tvd_limiter(scheme,cfl,fuu,fu,fd)
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
   REALTYPE,intent(in) :: cfl,fuu,fu,fd
!
! !LOCAL VARIABLES:
   REALTYPE           :: ratio,x,Phi
   REALTYPE,parameter :: one6th=_ONE_/6
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC

   if (abs(fd-fu) .gt. 1.d-10) then
      ratio = (fu-fuu)/(fd-fu)   ! slope ratio
   else
      ratio = (fu-fuu)*1.d10
   end if

   select case (scheme)
      case ((P2),(P2_PDM))
         x = one6th*(_ONE_-_TWO_*cfl)
         Phi = (_HALF_+x) + (_HALF_-x)*ratio
         if (scheme.eq.P2) then
            adv_tvd_limiter = Phi
         else
            adv_tvd_limiter = max(_ZERO_,min(Phi,_TWO_/(_ONE_-cfl),_TWO_*ratio/(cfl+1.d-10)))
         end if
      case (SUPERBEE)
         adv_tvd_limiter = max(_ZERO_,min(_ONE_,_TWO_*ratio),min(ratio,_TWO_))
      case (MUSCL)
         adv_tvd_limiter = max(_ZERO_,min(_TWO_,_TWO_*ratio,_HALF_*(_ONE_+ratio)))
      case default
         stop 'adv_tvd_limiter: invalid scheme'
   end select

   return
   end function adv_tvd_limiter
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
