!$Id: short_wave_radiation.F90,v 1.2 2003-04-23 12:05:50 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Short wave radiation.
!
! !INTERFACE:
   REALTYPE function short_wave_radiation(yday,hour,lat,lon,cc)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Short wave radiation is calculated based on the following input
!  parameters - year day, hour of day, latitude and longitude and cloud cover.
!  The albedo monthly values are from Payne (1972) as means of the values
!  at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
!  band of the Mediterranean Sea ) :
!  The radiation is returned as $W/m^2$.
!
! !SEE ALSO:
!  meteo.F90
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: yday
   REALTYPE, intent(inout)             :: hour
   REALTYPE, intent(in)                :: lat,lon,cc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding and Hans Burchard
!
!  $Log: short_wave_radiation.F90,v $
!  Revision 1.2  2003-04-23 12:05:50  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:39  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/07/26 14:35:18  bbh
!  initial import into CVS
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
   REALTYPE, parameter       :: deg2rad=pi/180.,rad2deg=180./pi

   REALTYPE, parameter       :: solar=1350.
   REALTYPE, parameter       :: eclips=23.439*deg2rad
   REALTYPE, parameter       :: tau=0.7
   REALTYPE, parameter       :: aozone=0.09

   REALTYPE, parameter       :: yrdays(2)= (/365.,366./)

   REALTYPE, parameter       :: alb1(20) = &
                    (/.719,.656,.603,.480,.385,.300,.250,.193,.164, &
                      .131,.103,.084,.071,.061,.054,.039,.036,.032,.031,.030 /)

   REALTYPE, parameter  :: za(20) = (/90.,88.,86.,84.,82.,80.,78.,76.,74.,70., &
                                      66.,62.,58.,54.,50.,40.,30.,20.,10.,0.0 /)
   REALTYPE                 :: dza(19)
   data            dza/8*2.0, 6*4.0, 5*10.0/

! !LOCAL VARIABLES:
   REALTYPE                  :: alat,alon,eqnx
   REALTYPE                  :: th0,th02,th03,sundec
   REALTYPE                  :: thsun,coszen,zen,dzen,sunbet
   REALTYPE                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   REALTYPE                  :: albedo
   integer                   :: jab
!
! !TO DO:
!
! !BUGS:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   th0 = 2.*pi*yday/yrdays(1)
   th02 = 2.*th0
   th03 = 3.*th0

!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)           &
           - 0.006758*cos(th02) + 0.000907*sin(th02)                   &
           - 0.002697*cos(th03) + 0.001480*sin(th03)

   alon = deg2rad*lon
   alat = deg2rad*lat

!  sun hour angle :
   thsun = (hour-12.)*15.*deg2rad + alon

!  cosine of the solar zenith angle :
   coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1./coszen)
   end if
   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

   eqnx = (yday-81.)/yrdays(1)*2.*pi

!  sin of the solar noon altitude in radians :
   sunbet=sin(alat)*sin(eclips*sin(eqnx))+cos(alat)*cos(eclips*sin(eqnx))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  calculates the albedo as a function of the solar zenith angle :
!  (after Payne jas 1972)
!  solar zenith angle in degrees :
   zen=(180./pi)*acos(coszen)
   if(zen .ge. 74.)then
      jab=.5*(90.-zen)+1.
   else if (zen .ge. 50.) then
      jab=.23*(74.-zen)+9.
   else
      jab=.10*(50.-zen)+15.
   endif

   dzen=(za(jab)-zen)/dza(jab)
   albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8

!HB I do not know where this "if(cc .lt. 0.3) then" comes from. It is not
!HB included in Rosati,Miyakoda 1988 ; eq. 3.8 ...
!HB   if(cc .lt. 0.3) then
!HB      short_wave_radiation  = qtot
!HB   else
      short_wave_radiation  = qtot*(1-0.62*cc + .0019*sunbet)*(1.-albedo)
!HB   endif

   return
   end function short_wave_radiation
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding & Hans Burchard
!-----------------------------------------------------------------------
