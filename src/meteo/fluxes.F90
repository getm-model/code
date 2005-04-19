!$Id: fluxes.F90,v 1.9 2005-04-19 12:21:23 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momentum fluxes.
!
! !INTERFACE:
   subroutine fluxes(u10,v10,airt,tcc,sst,hf,taux,tauy)
!
! !DESCRIPTION:
!  The sum of the latent and sensible heat fluxes + longwave
!  back-radiation is calculated and returned in \emph{hf} [$W/m^2$]. Also the
!  sea surface stresses are calculated and returned in \emph{taux} and
!  \emph{tauy} [$N/m^2$], repsectively. The wind velocities are following the
!  meteorological convention (from where) and are in $m/s$. The
!  temperatures \emph{airt} and \emph{sst} can be in Kelvin or Celcius -
!  if they are $>$ 100 - Kelvin is assumed. \emph{tcc} - the total cloud 
!  cover is specified as fraction between 0 and 1.
!
! !SEE ALSO:
!  meteo.F90, exchange_coefficients.F90
!
! !USES:
#define EA_ZERO
   use meteo, only: cpa,emiss,bolz,KELVIN
#ifdef EA_ZERO
   use meteo, only: w,L,rho_air,qs,qa
#else
   use meteo, only: w,L,rho_air,qs,qa,ea
#endif
   use meteo, only: cd_mom,cd_heat,cd_latent
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: u10,v10,airt,tcc,sst
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: hf,taux,tauy
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding and Hans Burchard
!
!  $Log: fluxes.F90,v $
!  Revision 1.9  2005-04-19 12:21:23  kbk
!  simulate old compiler bug by -DEA_ZERO
!
!  Revision 1.8  2005/01/13 09:49:37  kbk
!  wet bulb works, es is global, cleaning - Stips
!
!  Revision 1.7  2003/12/16 17:16:13  kbk
!  double declaration of cd_heat,cd_mom
!
!  Revision 1.6  2003/10/01 12:10:05  kbk
!  hum_method=1 (specific humidity) now works correctly
!
!  Revision 1.5  2003/07/01 16:38:34  kbk
!  cleaned code - new methods
!
!  Revision 1.4  2003/06/17 14:53:28  kbk
!  default meteo variables names comply with Adolf Stips suggestion + southpole(3)
!
!  Revision 1.3  2003/04/23 12:05:50  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/03/17 15:04:15  gotm
!  Fixed Kondo coefficients - -DWRONG_KONDO can be used
!
!  Revision 1.1.1.1  2002/05/02 14:01:39  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
!
!
! !DEFINED PARAMETERS:
   integer, parameter   :: clark=1      ! Clark et. al, 1974
   integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
!
! !LOCAL VARIABLES:
#ifdef EA_ZERO
   REALTYPE                  :: ea=_ZERO_
#endif
   REALTYPE                  :: tmp
   REALTYPE                  :: qe,qh,qb
   REALTYPE                  :: ta,ta_k,tw,tw_k
   integer                   :: back_radiation_method=clark
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+KELVIN
   else
      tw  = sst-KELVIN
      tw_k= sst
   end if

   if (airt .lt. 100.) then
      ta = airt
      ta_k  = airt + KELVIN
   else
      ta = airt - KELVIN
      ta_k = airt
   end if

   qh=cd_heat*cpa*rho_air*w*(tw-ta)            ! sensible
   qe=cd_latent*L*rho_air*w*(qs-qa)            ! latent

   select case(back_radiation_method)          ! back radiation
      case(clark)
!        AS - unit of ea is Pascal, must hPa
         qb=(1.0-.8*tcc*tcc)                                    &
            *emiss*bolz*(tw_k**4)*(0.39-0.05*sqrt(0.01*ea))     &
            +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case(hastenrath) ! qa in g(water)/kg(wet air)
         qb=(1.0-.8*tcc*tcc)                                    &
            *emiss*bolz*(tw_k**4)*(0.39-0.056*sqrt(1000.0*qa))  &
            +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case default
   end select

   hf = -(qe+qh+qb)

   tmp   = cd_mom*rho_air*w
   taux  = tmp*u10
   tauy  = tmp*v10

   return
   end subroutine fluxes
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
