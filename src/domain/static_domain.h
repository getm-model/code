#include "dimensions.h"

!  coordinate axes - grid-type = 1 or 2
   REALTYPE,target                     :: xcord(imin-HALO:imax+HALO)
   REALTYPE,target                     :: ycord(jmin-HALO:jmax+HALO)

!  pseudo coordinate axes - grid-type = 3 or 4
   REALTYPE,target                     :: xxcord(imin-HALO-1:imax+HALO)
   REALTYPE,target                     :: yxcord(jmin-HALO-1:jmax+HALO)

!  mask
   REALTYPE                            :: mask(E2DFIELD)
   integer,target                      :: az(E2DFIELD)
   integer,target                      :: au(E2DFIELD)
   integer,target                      :: av(E2DFIELD)
   integer                             :: ax(E2DFIELD)

!  bathymetry
   REALTYPE                            :: H(E2DFIELD)  = -10.
   REALTYPE                            :: HU(E2DFIELD) = -10.
   REALTYPE                            :: HV(E2DFIELD) = -10.
   REALTYPE                            :: dry_z(E2DFIELD)
   REALTYPE                            :: dry_u(E2DFIELD)
   REALTYPE                            :: dry_v(E2DFIELD)

!  coriolis terms
   REALTYPE                            :: cor(E2DFIELD)
   REALTYPE                            :: coru(E2DFIELD)
   REALTYPE                            :: corv(E2DFIELD)

!  lat/lon
   REALTYPE,target                     :: lonc(E2DFIELD) = -999.
   REALTYPE,target                     :: latc(E2DFIELD) = -999.
   REALTYPE,target                     :: lonx(E2DXFIELD) = -999.
   REALTYPE,target                     :: latx(E2DXFIELD) = -999.
   REALTYPE                            :: lonu(E2DFIELD) = -999.
   REALTYPE                            :: latu(E2DFIELD) = -999.
   REALTYPE                            :: lonv(E2DFIELD) = -999.
   REALTYPE                            :: latv(E2DFIELD) = -999.

!  grid convergence
!KB   REALTYPE                            :: angle(E2DFIELD)
   REALTYPE                            :: convc(E2DFIELD) = _ZERO_
   REALTYPE                            :: convx(E2DXFIELD) = _ZERO_

!  grid points
   REALTYPE,target                     :: xx(E2DXFIELD)
   REALTYPE,target                     :: yx(E2DXFIELD)
   REALTYPE,target                     :: xc(E2DFIELD)
   REALTYPE,target                     :: yc(E2DFIELD)
   REALTYPE                            :: xu(E2DFIELD)
   REALTYPE                            :: yu(E2DFIELD)
   REALTYPE                            :: xv(E2DFIELD)
   REALTYPE                            :: yv(E2DFIELD)

!  metric parameters
   REALTYPE                            :: dx=-_ONE_,dy=-_ONE_,x0,y0,ard1
   REALTYPE                            :: dlon=-_ONE_,dlat=-_ONE_,lon0,lat0
   REALTYPE                            :: dxdyc(E2DFIELD)
   REALTYPE                            :: dydxc(E2DFIELD)
   REALTYPE,target                     :: dxc(E2DFIELD) = -999.
   REALTYPE,target                     :: dxu(E2DFIELD) = -999.
   REALTYPE,target                     :: dxv(E2DFIELD) = -999.
   REALTYPE,target                     :: dxx(E2DFIELD) = -999.
   REALTYPE,target                     :: dyc(E2DFIELD) = -999.
   REALTYPE,target                     :: dyu(E2DFIELD) = -999.
   REALTYPE,target                     :: dyv(E2DFIELD) = -999.
   REALTYPE,target                     :: dyx(E2DFIELD) = -999.
   REALTYPE,target                     :: arcd1(E2DFIELD)
   REALTYPE,target                     :: arud1(E2DFIELD)
   REALTYPE,target                     :: arvd1(E2DFIELD)

!  bottom roughness
   REALTYPE,dimension(E2DFIELD)        :: z0,zub0,zvb0
