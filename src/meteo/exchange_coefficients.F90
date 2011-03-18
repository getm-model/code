#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Air/sea exchange coefficients
!
! !INTERFACE:
   subroutine exchange_coefficients(u10,v10,airt,airp,sst,hum,hum_method)
!
! !DESCRIPTION:
!  Various variables for calculating meteorological forcing is calculated
!  here. The input variables are given in SI units ($m/s$, $^oC | Kelvin$,
!  and $Pa$). The scheme used to get the required variables is as follows:
!  \begin{itemize}
!     \item We have sst, wind, air temperature
!     \item We calculate the saturation vapor pressure (svp(T) based on 
!           the sea surface temperature.
!     \item We calculate the specific humidity valid for the sst - using
!           the just calculated the saturation vapor pressure - unit must
!           be [kg/kg]
!     \item We have some measure of the humidity of the air - either
!           specific humidity at air temperature, relative humidity or
!           wet bulb tmperature. This needs to be converted to specific 
!           humidity at 2m.
!     \item In the case of precipitation, we compute the sensible heatflux
!           due to the additional water, assuming that the rain has the same
!           temperature as the air.
!  \end{itemize}
!  
!  The following formulaes are used (for the saturation vapor pressure
!  a large number of different formulaes exists):
!  \begin{itemize}
!     \item Latent heat: $L = (2.5-0.00234T)10^6$
!     \item Saturation vapor pressure: $e_s = a_0 + a_1 T^1 + a_2 T^2 + a_3 T^3
!                                    + a_4 T^4 + a_5 T^5 + a_6 T^6 + a_7 T^7$
!     \item Specific humidity: $q_s = const06 e_s / (airp-0.377 e_s) $
!     \item Actual humidity: $q_a = 0.01  rh q_s $
!     \item Virtual temperature: $T_v = T_k ( 1 + q_a/const06 )/( 1 + q_a)$
!     \item Air density: $\rho_a = airp/(287.05 T_v)$
!     \item Air stability: $S_0 = 0.25(t_w - t_a)/W and s = S_0 |S_0|
!                                                         / (|S_0| +0.01)$
!  \end{itemize}
!  The variables $L, e_s, q_s, e_a$ and $q_a$ are later used for calculating
!  the latent and sensible heat fluxes.
!  The exchange coefficients $C_{D_{mom}}$ and $C_{D_{heat}}$ are functions of
!  emperical parameters - which a actual wind speed dependend - and the
!  air stability.
!
! !SEE ALSO:
!  meteo.F90, fluxes.F90
!
! !USES:
   use meteo, only: cpa,KELVIN
   use meteo, only: L,rho_air,w,qs,qa,ea,es
   use meteo, only: cd_mom,cd_heat,cd_latent
   use meteo, only: fwf_method,cpw,cd_precip
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: u10,v10,airt,airp,sst,hum
   integer, intent(in)                 :: hum_method
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: a1=6.107799961
   REALTYPE, parameter       :: a2=4.436518521e-1
   REALTYPE, parameter       :: a3=1.428945805e-2
   REALTYPE, parameter       :: a4=2.650648471e-4
   REALTYPE, parameter       :: a5=3.031240396e-6
   REALTYPE, parameter       :: a6=2.034080948e-8
   REALTYPE, parameter       :: a7=6.136820929e-11
   REALTYPE, parameter       :: const06=0.62198
   REALTYPE, parameter       :: eps=1.0e-12
   REALTYPE, parameter       :: rgas=287.1  !J/kg/K gas const dry air
!
! !LOCAL VARIABLES:
   REALTYPE                  :: tvirt,s,s0
   REALTYPE                  :: ae_d,be_d,ce_d,pe_d
   REALTYPE                  :: ae_h,be_h,ce_h,pe_h
   REALTYPE                  :: ae_e,be_e,ce_e,pe_e
   REALTYPE                  :: cee_d,cee_h,cee_e
   REALTYPE                  :: x
   REALTYPE                  :: ta,ta_k,tw,tw_k

   REALTYPE                  :: rh
   REALTYPE                  :: twet,twet_k
   REALTYPE                  :: dew,dew_k
   REALTYPE                  :: x1,x2,x3
!EOP
!-----------------------------------------------------------------------
!BOC

! BJB-TODO: Change paramenters+constants in code to double precision

!  water temperature
   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+KELVIN
   else
      tw  = sst-KELVIN
      tw_k= sst
   end if

!  air temperature
   if (airt .lt. 100.) then
      ta   = airt
      ta_k = airt+KELVIN
   else
      ta   = airt-KELVIN
      ta_k = airt
   end if

!  windspeed
   w = sqrt(u10*u10+v10*v10)
!  latent heat of vaporisation
   L = (2.5-0.00234*tw)*1.e6

!  saturation vapor pressure - using SST as temperature
!  various formulations are available
!  I've tested in a spread sheet - they give almost the same
!  result
!  for hum_method=2 the calculation of es can be substituted
!  by any of the formlaes below.
#if 0
!  http://www.cdc.noaa.gov/coads/software/other/profs_short
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal

!  http://www.cdc.noaa.gov/coads/software/other/profs_short
   es = 100.*const06*exp(17.67*tw/(tw+243.5))

!  http://www.usatoday.com/weather/whumcalc.htm
   es = 6.11*10.0**(7.5*tw/(237.7+tw))

!  From the HAMSOM model
   es  = 611.21*exp((18.729 - (min(tw_k,300.)-273.15)/227.3)*       &
         (min(tw_k,300.)-273.15)/(max(tw_k,200.)-273.15+257.87))
#endif

#ifdef HAMSOM_SVP
    x  = (18.729 - (min(tw_k,300.)-273.15)/227.3)*       &
         (min(tw_k,300.)-273.15)/(max(tw_k,240.)-273.15+257.87)
    es = 611.2 * exp(x)
#else
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal
#endif
#ifndef OLD_WRONG_FLUXES
!  correction for seawater, following Kraus 1972
!  correcting for salt water assuming 98% RH
   es=0.98 * es
#endif
!  saturation specific humidity
   qs = const06*es/(airp-0.377*es)

!  see ../ncdf/ncdf_meteo.F90 for defined constants
   select case (hum_method)
   case (1)
!     specific humidity in kg/kg is given
      qa = hum
!     actual water vapor pressure in Pascal
      ea = qa *airp/(const06+0.378*qa)
   case (2)
!     relative humidity in % given
      rh = 0.01 * hum
!     saturation vapor pressure at that air temperature
#ifdef HAMSOM_SVP
      x1 = (18.729 - (min(ta_k,300.*_ONE_)-273.15)/227.3)
      x2 = (min(ta_k,300.*_ONE_)-273.15)
      x3 = (max(ta_k,200.*_ONE_)-273.15+257.87)
      ea = 611.21*exp(x1*x2/x3)
#else
      ea = a1 +ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*(a6+ta*a7)))))
      ea = ea * 100.0 ! Conversion millibar --> Pascal
#endif
!     get actual vapor pressure
      ea = rh * ea
!     convert to specific humidity
      qa = const06*ea/(airp-0.377*ea)
   case (3)
      ! Piece of code taken from HAMSOM for calculating relative
      ! humidity from dew point temperature and dry air temperature.
      ! It must be sure that hum is dew point temperature in Kelvin 
      ! in the next line ...
!       use dew in degC
      if (hum .lt. 100.) then
         dew = hum
         dew_k = hum+KELVIN
      else
         dew = hum-KELVIN
         dew_k = hum
      end if
#ifdef HAMSOM_SVP
      x1 = (18.729 - (min(dew_k,300.*_ONE_)-273.15)/227.3)
      x2 = (min(dew_k,300.*_ONE_)-273.15)
      x3 = (max(dew_k,200.*_ONE_)-273.15+257.87)
      ea = 611.21*exp(x1*x2/x3)
#ifdef OLD_WRONG_FLUXES
      x1 = (18.729 - (min(ta_k,300.*_ONE_)-273.15)/227.3)
      x2 = (min(ta_k,300.*_ONE_)-273.15)
      x3 = (max(ta_k,200.*_ONE_)-273.15+257.87)
      es = 611.21*exp(x1*x2/x3)
#endif
#else
      ea = a1 +dew*(a2+dew*(a3+dew*(a4+dew*(a5+dew*(a6+dew*a7)))))
      ea = ea * 100.0 ! Conversion millibar --> Pascal
#endif

#ifdef OLD_WRONG_FLUXES
      rh = 100.*ea/es
      qa = 0.01*rh*qs
#else
!     AS convert actual vapor pressure to specific humidity
      qa = const06*ea/(airp-0.377*ea)
#endif
   case (4)
!     wet bulb temperature given
!     Calculate the SVP at wet bulb temp then
!     use the psychrometer formula to get vapour pressure
!     See Smithsonian Met tables 6th Edition pg 366 eqn 3
!     Make sure this is in degC 
      if (hum .lt. 100 ) then
         twet=hum 
         twet_k=hum+KELVIN
      else
         twet=hum-KELVIN
         twet_k=hum
      end if
!     saturation vapor pressure at wet bulb temperature
#ifdef HAMSOM_SVP
      x  = (18.729 - (min(twet_k,300.)-273.15)/227.3)*       &
           (min(twet_k,300.)-273.15)/(max(twet_k,240.)-273.15+257.87)
      ea = 611.2 * exp(x)
#else
      ea = a1 +twet*(a2+twet*(a3+twet*(a4+twet*(a5+twet*(a6+twet*a7)))))
      ea = ea * 100.0 ! Conversion millibar --> Pascal
#endif
!     actual vapor pressure
      ea = ea - 6.6e-4*(1+1.15e-3*twet)*airp*(ta-twet)
!     specific humidity in kg/kg
      qa = const06*ea/(airp-0.377*ea)
   case default
      FATAL 'not a valid hum_method'
      stop 'exchange_coefficients()'
   end select

   tvirt = ta_k*(1+qa/const06)/(1+qa)
   rho_air = airp/(287.05*Tvirt)

!  Transfer coefficient for heat and momentum
   if (w .lt. 2.2) then
      ae_d=0.0;   be_d=1.08;   ce_d=0.0;      pe_d=-0.15;
      ae_h=0.0;   be_h=1.185;  ce_h=0.0;      pe_h=-0.157;
      ae_e=0.0;   be_e=1.23;   ce_e=0.0;      pe_e=-0.16;
   else if (w .lt. 5.0) then
      ae_d=0.771; be_d=0.0858; ce_d=0.0;      pe_d=1.0;
      ae_h=0.927; be_h=0.0546; ce_h=0.0;      pe_h=1.0;
      ae_e=0.969; be_e=0.052;  ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 8.0) then
      ae_d=0.867; be_d=0.0667; ce_d=0.0;      pe_d=1.0;
      ae_h=1.15;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_e=1.18;  be_e=0.01;   ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 25.0) then
      ae_d=1.2;   be_d=0.025;  ce_d=0.0;      pe_d=1.0;
      ae_h=1.17;  be_h=0.0075; ce_h=-0.00045; pe_h=1.0
      ae_e=1.196; be_e=0.008;  ce_e=-0.0004;  pe_e=1.0
   else
      ae_d=0.0;   be_d=0.073;  ce_d=0.0;      pe_d=1.0;
      ae_h=1.652; be_h=-0.017; ce_h=0.0;      pe_h=1;
      ae_e=1.68;  be_e=-0.016; ce_e=0.0;      pe_e=1;
   end if

   cee_d = (ae_d+be_d*exp(pe_d*log(w+eps)))*1.0e-3
   cee_h = (ae_h+be_h*exp(pe_h*log(w+eps))+ce_h*(w-8.0)**2)*1.0e-3
   cee_e = (ae_e+be_e*exp(pe_e*log(w+eps))+ce_e*(w-8.0)**2)*1.0e-3

!  Stability - assume 10 meter wind
   s0=0.25*(tw-ta)/(w+eps)**2
   s=s0*abs(s0)/(abs(s0)+0.01)

   if(s .lt. -3.3) then
      cd_mom = 0.0 ; cd_heat = 0.0; cd_latent = 0.0
   else if(s .lt. 0.0) then
      x = 0.1+0.03*s+0.9*exp(4.8*s)
      cd_mom =x*cee_d
      cd_heat=x*cee_h
      cd_latent=x*cee_e
   else
      cd_mom =cee_d*(1.0+0.47*sqrt(s))
      cd_heat=cee_h*(1.0+0.63*sqrt(s))
      cd_latent=cee_e*(1.0+0.63*sqrt(s))
   end if

   if (fwf_method .ge. 1) then
      x1 = 2.11e-5*(ta_k/KELVIN)**1.94
      x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rho_air*cpa)
      x3 = qa * L /(rgas * ta_K * ta_K)
      cd_precip = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
      cd_precip = cd_precip*cpw*((tw-ta) + (qs-qa)*L/cpa)
   else
      cd_precip = _ZERO_
   end if

   return
   end subroutine exchange_coefficients
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
