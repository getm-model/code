#include "dimensions.h"

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
   REALTYPE                            :: lonc(E2DFIELD)
   REALTYPE                            :: latc(E2DFIELD)
   REALTYPE                            :: lonx(E2DFIELD)
   REALTYPE                            :: latx(E2DFIELD)
   REALTYPE                            :: lonu(E2DFIELD)
   REALTYPE                            :: latu(E2DFIELD)
   REALTYPE                            :: lonv(E2DFIELD)
   REALTYPE                            :: latv(E2DFIELD)

!  grid convergence
   REALTYPE                            :: angle(E2DFIELD)
   REALTYPE                            :: convc(E2DFIELD)
   REALTYPE                            :: convx(E2DFIELD)

!  grid points
   REALTYPE                            :: xx(E2DFIELD)
   REALTYPE                            :: yx(E2DFIELD)
   REALTYPE                            :: xc(E2DFIELD)
   REALTYPE                            :: yc(E2DFIELD)
   REALTYPE                            :: xu(E2DFIELD)
   REALTYPE                            :: yu(E2DFIELD)
   REALTYPE                            :: xv(E2DFIELD)
   REALTYPE                            :: yv(E2DFIELD)

!  metric parameters
   REALTYPE                            :: dx,dy,x0,y0,ard1
   REALTYPE                            :: dlon,dlat,lon0,lat0
   REALTYPE                            :: dxdyc(E2DFIELD)
   REALTYPE                            :: dydxc(E2DFIELD)
   REALTYPE                            :: dxc(E2DFIELD)
   REALTYPE                            :: dxu(E2DFIELD)
   REALTYPE                            :: dxv(E2DFIELD)
   REALTYPE                            :: dxx(E2DFIELD)
   REALTYPE                            :: dyc(E2DFIELD)
   REALTYPE                            :: dyu(E2DFIELD)
   REALTYPE                            :: dyv(E2DFIELD)
   REALTYPE                            :: dyx(E2DFIELD)
   REALTYPE                            :: arcd1(E2DFIELD)
   REALTYPE                            :: arud1(E2DFIELD)
   REALTYPE                            :: arvd1(E2DFIELD)

