#include "dimensions.h"

!  coordinate axes - grid-type = 1 or 2
   REALTYPE                             :: xcord(imin-HALO:imax+HALO)
   REALTYPE                             :: ycord(jmin-HALO:jmax+HALO)

!  pseudo coordinate axes - grid-type = 3 or 4
   REALTYPE                             :: xxcord(imin-HALO-1:imax+HALO)
   REALTYPE                             :: yxcord(jmin-HALO-1:jmax+HALO)

!  mask
   integer                             :: az(E2DFIELD)
   integer                             :: au(E2DFIELD)
   integer                             :: av(E2DFIELD)
   integer                             :: ax(E2DFIELD)

!  bathymetry
   REALTYPE                            :: H(E2DFIELD)
   REALTYPE                            :: HU(E2DFIELD)
   REALTYPE                            :: HV(E2DFIELD)
   REALTYPE                            :: dry_z(E2DFIELD)
   REALTYPE                            :: dry_u(E2DFIELD)
   REALTYPE                            :: dry_v(E2DFIELD)

!  coriolis terms
   REALTYPE                            :: cor(E2DFIELD)
   REALTYPE                            :: coru(E2DFIELD)
   REALTYPE                            :: corv(E2DFIELD)

!  lat/lon
   REALTYPE                            :: lonc(E2DFIELD) = -999.
   REALTYPE                            :: latc(E2DFIELD) = -999.
   REALTYPE                            :: lonx(E2DXFIELD) = -999.
   REALTYPE                            :: latx(E2DXFIELD) = -999.
   REALTYPE                            :: lonu(E2DFIELD) = -999.
   REALTYPE                            :: latu(E2DFIELD) = -999.
   REALTYPE                            :: lonv(E2DFIELD) = -999.
   REALTYPE                            :: latv(E2DFIELD) = -999.

!  grid convergence
   REALTYPE                            :: angle(E2DFIELD)
   REALTYPE                            :: convc(E2DFIELD) = -999.
   REALTYPE                            :: convx(E2DXFIELD) = -999.

!  grid points
   REALTYPE                            :: xx(E2DXFIELD)
   REALTYPE                            :: yx(E2DXFIELD)
   REALTYPE                            :: xc(E2DFIELD)
   REALTYPE                            :: yc(E2DFIELD)
   REALTYPE                            :: xu(E2DFIELD)
   REALTYPE                            :: yu(E2DFIELD)
   REALTYPE                            :: xv(E2DFIELD)
   REALTYPE                            :: yv(E2DFIELD)

!  metric parameters
   REALTYPE                            :: dx=-_ONE_,dy=-_ONE_,x0,y0,ard1
   REALTYPE                            :: dlon=-_ONE_,dlat=-_ONE_,lon0,lat0
   REALTYPE                            :: dxdyc(E2DFIELD)
   REALTYPE                            :: dydxc(E2DFIELD)
   REALTYPE                            :: dxc(E2DFIELD) = -999.
   REALTYPE                            :: dxu(E2DFIELD) = -999.
   REALTYPE                            :: dxv(E2DFIELD) = -999.
   REALTYPE                            :: dxx(E2DFIELD) = -999.
   REALTYPE                            :: dyc(E2DFIELD) = -999.
   REALTYPE                            :: dyu(E2DFIELD) = -999.
   REALTYPE                            :: dyv(E2DFIELD) = -999.
   REALTYPE                            :: dyx(E2DFIELD) = -999.
   REALTYPE                            :: arcd1(E2DFIELD)
   REALTYPE                            :: arud1(E2DFIELD)
   REALTYPE                            :: arvd1(E2DFIELD)

!  bottom roughness
   REALTYPE                            :: z0(E2DFIELD)
