   integer		:: iextr=-1, jextr=-1
   integer		:: imin=-1,imax=-1,jmin=-1,jmax=-1
   integer		:: iimin=-1,iimax=-1,jjmin=-1,jjmax=-1
   integer		:: kmax=1
   REALTYPE,dimension(:,:),allocatable	:: H,HU,HV
   REALTYPE,dimension(:,:),allocatable	:: lonmap,latmap,conv
   REALTYPE,dimension(:,:),allocatable	:: dry_z,dry_u,dry_v
   REALTYPE,dimension(:,:),allocatable	:: cor,coru,corv
#if ! ( defined(CURVILINEAR) || defined(SPHERICAL) )
   REALTYPE,dimension(:),allocatable	:: xc,yc
   REALTYPE				:: dx,dy,ard1
#else 
   REALTYPE,dimension(:,:),allocatable	:: xc,yc,xu,yu,xv,yv,xx,yx
   REALTYPE,dimension(:,:),allocatable	:: lonx,latx,lonc,latc
   REALTYPE,dimension(:,:),allocatable	:: lonu,latu,lonv,latv
   REALTYPE,dimension(:,:),allocatable	:: dxdyc,dydxc,angle
   REALTYPE,dimension(:,:),allocatable	:: dxc,dxu,dxv,dxx
   REALTYPE,dimension(:,:),allocatable	:: dyc,dyu,dyv,dyx
   REALTYPE,dimension(:,:),allocatable	:: arcd1,arud1,arvd1
#endif
   integer, dimension(:,:),allocatable	:: az,au,av,ax

