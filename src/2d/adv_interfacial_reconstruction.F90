#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  REALTYPE function adv_interfacial_reconstruction -
!
! !INTERFACE:
   REALTYPE function adv_interfacial_reconstruction(scheme,cfl,fuu,fu,fd)
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
   REALTYPE           :: ratio,limiter,x,deltaf,deltafu
   REALTYPE,parameter :: one6th=_ONE_/6
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC

   deltaf  = fd - fu
   deltafu = fu - fuu

   if (scheme .eq. P2) then
      x = one6th*(_ONE_-_TWO_*cfl)
      adv_interfacial_reconstruction = fu + _HALF_*(_ONE_-cfl)*((_HALF_+x)*deltaf + (_HALF_-x)*deltafu)
      return
   end if

!  Upstream
   adv_interfacial_reconstruction = fu

   if (deltaf*deltafu .gt. _ZERO_) then

      ratio = deltafu / deltaf   ! slope ratio

      select case (scheme)
         case (P2_PDM)
            x = one6th*(_ONE_-_TWO_*cfl)
            limiter = (_HALF_+x) + (_HALF_-x)*ratio
!            limiter = max(_ZERO_,min(_TWO_*ratio/(cfl+1.d-10),limiter,_TWO_/(_ONE_-cfl)))
            limiter = min(_TWO_*ratio/(cfl+1.d-10),limiter,_TWO_/(_ONE_-cfl))
         case (SUPERBEE)
!            limiter = max(_ZERO_,min(_TWO_*ratio,_ONE_),min(ratio,_TWO_))
            limiter = max(min(_TWO_*ratio,_ONE_),min(ratio,_TWO_))
         case (MUSCL)
!            limiter = max(_ZERO_,min(_TWO_*ratio,_HALF_*(_ONE_+ratio),_TWO_))
            limiter = min(_TWO_*ratio,_HALF_*(_ONE_+ratio),_TWO_)
         case default
            return
      end select

      adv_interfacial_reconstruction = fu + _HALF_*limiter*(_ONE_-cfl)*deltaf

   end if

   return
   end function adv_interfacial_reconstruction
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
