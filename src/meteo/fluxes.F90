!$Id: fluxes.F90,v 1.1 2002-05-02 14:01:39 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momentum fluxes.
!
! !INTERFACE:
   subroutine fluxes(u10,v10,airt,cc,sst,hf,taux,tauy)
!
! !DESCRIPTION:
!  The sum of the latent and sensible heat fluxes + longwave
!  back-radiation is calculated and returned in \emph{hf} [$W/m^2$]. Also the
!  sea surface stresses are calculated and returned in \emph{taux} and
!  \emph{tauy} [$N/m^2$], repsectively. The wind velocities are following the
!  meteorological convention (from where) and are in $m/s$. The
!  temperatures \emph{airt} and \emph{sst} can be in Kelvin or Celcius -
!  if they are $>$ 100 - Kelvin is assumed. \emph{cc} - the cloud cover -
!  is specified as fraction between 0 and 1.
!
! !SEE ALSO:
!  meteo.F90, exchange_coefficients.F90
!
! !USES:
   use meteo, only: cpa,emiss,bolz,KELVIN
   use meteo, only: w,L,rho_air,qs,qa,ea,cd_heat,cd_mom
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)		:: u10,v10,airt,cc,sst
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: hf,taux,tauy
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding and Hans Burchard
!
!  $Log: fluxes.F90,v $
!  Revision 1.1  2002-05-02 14:01:39  gotm
!  Initial revision
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
   REALTYPE		:: tmp
   REALTYPE		:: qe,qh,qb
   REALTYPE		:: ta,tw,tw_k
   integer		:: back_radiation_method=clark
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

   if (airt .gt. 100.) then
      ta  = airt - KELVIN
   else
      ta = airt
   end if

   qe=cd_heat*L*rho_air*w*(qs-qa)			! latent
   qh=cd_heat*cpa*rho_air*w*(tw-ta)			! sensible

   select case(back_radiation_method)			! back radiation
      case(clark)
         qb=(1.0-.8*cc*cc)				&
            *emiss*bolz*(tw_k**4)*(0.39-0.05*sqrt(ea/100.0))	&
            +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case(hastenrath) ! qa in g(water)/kg(wet air)
         qb=(1.0-.8*cc*cc)				&
            *emiss*bolz*(tw_k**4)*(0.39-0.056*sqrt(1000*qa))		&
            +4.0*emiss*bolz*(tw_k**3)*(tw-ta)
      case default
   end select

#if 0
   hf = -qb
   hf = -qh
   hf = -qe
#endif
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
