   integer         :: iextr=-1, jextr=-1
   integer         :: imin=-1,imax=-1,jmin=-1,jmax=-1
   integer         :: kmax=1

!  mask
   integer, dimension(:,:), allocatable     :: az,au,av,ax

!  bathymetry
   REALTYPE, dimension(:,:), allocatable    :: H,HU,HV
   REALTYPE, dimension(:,:), allocatable    :: dry_z,dry_u,dry_v

!  coriolis terms
   REALTYPE, dimension(:,:), allocatable    :: cor,coru,corv

!  lat/lon
   REALTYPE, dimension(:,:), allocatable    :: lonc,latc
   REALTYPE, dimension(:,:), allocatable    :: lonx,latx
   REALTYPE, dimension(:,:), allocatable    :: lonu,latu
   REALTYPE, dimension(:,:), allocatable    :: lonv,latv

!  grid convergence
   REALTYPE, dimension(:,:), allocatable    :: angle,convc,convx

!  grid points
   REALTYPE, dimension(:,:), allocatable    :: xx,yx
   REALTYPE, dimension(:,:), allocatable    :: xc,yc
   REALTYPE, dimension(:,:), allocatable    :: xu,yu
   REALTYPE, dimension(:,:), allocatable    :: xv,yv

!  metric parameters
   REALTYPE                                 :: dx,dy,x0,y0,ard1
   REALTYPE                                 :: dlon,dlat,lon0,lat0
   REALTYPE, dimension(:,:), allocatable    :: dxdyc,dydxc
   REALTYPE, dimension(:,:), allocatable    :: dxc,dxu,dxv,dxx
   REALTYPE, dimension(:,:), allocatable    :: dyc,dyu,dyv,dyx
   REALTYPE, dimension(:,:), allocatable    :: arcd1,arud1,arvd1

!  bottom roughness
   REALTYPE, dimension(:,:), allocatable    :: z0
