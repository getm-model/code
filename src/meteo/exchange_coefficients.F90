!$Id: exchange_coefficients.F90,v 1.3 2003-04-23 12:05:50 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Air/sea exchange coefficients
!
! !INTERFACE:
   subroutine exchange_coefficients(u10,v10,airt,airp,sst,hum)
!
! !DESCRIPTION:
!  Various variables for calculating meteorological forcing is calculated
!  here. The input variables are gibe in SI units ($m/s$, $^oC | Kelvin$,
!  and $Pa$). The following formulaes are used:
!  \begin{itemize}
!     \item Latent heat: $L = (2.5-0.00234T)10^6$
!     \item Specific vapor pressure: $e_s = a_0 + a_1 T^1 + a_2 T^2 + a_3 T^3
!                                    + a_4 T^4 + a_5 T^5 + a_6 T^6 + a_7 T^7$
!     \item Specific humidity: $q_s = const06 e_s / (airp-0.377 e_s) $
!     \item Absolute humidity: $q_a = 0.01  rh q_s $
!     \item Absolute vapor pressure: $q_s = q_a airp / (airp + 0.377 q_a) $
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
#ifdef WRONG_KONDO
   use meteo, only: L,rho_air,w,qs,qa,cd_heat,cd_mom
#else
   use meteo, only: L,rho_air,w,qs,qa,cd_mom,cd_heat,cd_latent
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: u10,v10,airt,airp,sst,hum
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: exchange_coefficients.F90,v $
!  Revision 1.3  2003-04-23 12:05:50  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/03/17 15:04:15  gotm
!  Fixed Kondo coefficients - -DWRONG_KONDO can be used
!
!  Revision 1.1.1.1  2002/05/02 14:01:38  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
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
!
! !LOCAL VARIABLES:
   REALTYPE                  :: tvirt,s,s0
#ifdef WRONG_KONDO
   REALTYPE                  :: ae_h,be_h,ce_h,pe_h
   REALTYPE                  :: ae_m,be_m,ce_m,pe_m
   REALTYPE                  :: cee_heat,cee_mom
#else
   REALTYPE                  :: ae_d,be_d,ce_d,pe_d
   REALTYPE                  :: ae_h,be_h,ce_h,pe_h
   REALTYPE                  :: ae_e,be_e,ce_e,pe_e
   REALTYPE                  :: cee_d,cee_h,cee_e
#endif
   REALTYPE                  :: x,es,ea
   REALTYPE                  :: ta,ta_k,tw,tw_k

   REALTYPE                  :: twet,rh
#ifdef ECMWF_FRV
   REALTYPE                  :: d_tmp
#endif
!
!  !TO DO:
!
!  !BUGS:
!   How to calculate absolute humidity based on the input variables.
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
      ta   = airt
      ta_k = airt+KELVIN
   else
      ta   = airt-KELVIN
      ta_k = airt
   end if

   w = sqrt(u10*u10+v10*v10)
   L = (2.5-0.00234*tw)*1.e6
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal
   qs = const06*es/(airp-0.377*es)

!  FIXME
!HB   rh = 90.0
!  FIXME

#ifdef ECMWF_FRV
!     Piece of code taken from HAMSOM for calculating relative
!     humidity from dew point temperature and dry air temperature.
      d_tmp = hum
! It must be sure that hum is dew point temperature in Kelvin in the next
! line ...
      ea  = 611.21*exp((18.729 - (min(d_tmp,300.)-273.15)/227.3)*       &
              (min(d_tmp,300.)-273.15)/(max(d_tmp,200.)-273.15+257.87))
      es  = 611.21*exp((18.729 - (min(ta_k,300.)-273.15)/227.3)*       &
              (min(ta_k,300.)-273.15)/(max(ta_k,200.)-273.15+257.87))
      rh = ea/es * 100.
#endif

   if (rh .lt. 0.0) then
      ea = es - 67.*(ta-twet);
      x = (ta-twet)/(CONST06*L);
      ea = (es-cpa*airp*x)/(1+cpa*x);
      if(ea .lt. 0.0) ea = 0.0
      qa = CONST06*ea/(airp-0.377*ea);
   else
      qa = 0.01*rh*qs;
      ea = qa*airp/(const06 + 0.377*qa);
   end if

   tvirt = ta_k*(1+qa/const06)/(1+qa)
   rho_air = airp/(287.05*Tvirt)

!  Transfer coefficient for heat and momentum
#ifdef WRONG_KONDO
   if (w .lt. 2.2) then
      ae_h=0.0;   be_h=1.23;   ce_h=0.0;      pe_h=-0.16;
      ae_m=0.0;   be_m=1.185;  ce_m=0.0;      pe_m=-0.157;
   else if (w .lt. 5.0) then
      ae_h=0.969; be_h=0.0521; ce_h=0.0;      pe_h=1.0;
      ae_m=0.927; be_m=0.0546; ce_m=0.0;      pe_m=1.0;
   else if (w .lt. 8.0) then
      ae_h=1.18;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_m=1.15;  be_m=0.01;   ce_m=0.0;      pe_m=1.0;
   else if (w .lt. 25.0) then
      ae_h=1.196; be_h=0.008;  ce_h=-0.0004;  pe_h=1.0
      ae_m=1.17;  be_m=0.0075; ce_m=-0.00044; pe_m=1.0;
   else
      ae_h=1.68;  be_h=-0.016; ce_h=0;        pe_h=1;
      ae_m=1.652; be_m=-0.017; ce_m=0.0;      pe_m=1.0;
   end if

   cee_heat=(ae_h+be_h*exp(pe_h*log(w+eps))+ce_h*(w-8.0)**2)*1.0e-3
   cee_mom =(ae_m+be_m*exp(pe_m*log(w+eps)))*1.0e-3

!  Stability - assume 10 meter wind
   s0=0.25*(tw-ta)/(w+eps)**2
   s=s0*abs(s0)/(abs(s0)+0.01)

   if(s .lt. 0.) then
      x = 0.1+0.03*s+0.9*exp(4.8*s)
      cd_heat=x*cee_heat
      cd_mom =x*cee_mom
   else
      cd_heat=cee_heat*(1.0+0.63*sqrt(s))
      cd_mom =cee_mom *(1.0+0.47*sqrt(s))
   end if

#else
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
#endif

   return
   end subroutine exchange_coefficients
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
