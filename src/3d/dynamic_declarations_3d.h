! Remember to update this value if you add more 3D arrays.
#ifdef UV_TVD
  integer,parameter                    :: n3d_fields=33
#else
  integer,parameter                    :: n3d_fields=26
#endif

! Number of vertical layers in z,u,v columns
  integer, dimension(:,:), allocatable:: kmin,kumin,kvmin
  integer, dimension(:,:), allocatable:: kmin_pmz,kumin_pmz,kvmin_pmz

  REALTYPE, dimension(:,:,:), allocatable   :: uu,vv,ww
  REALTYPE, dimension(:,:,:), allocatable   :: ho,hn
  REALTYPE, dimension(:,:,:), allocatable   :: huo,hun
  REALTYPE, dimension(:,:,:), allocatable   :: hvo,hvn
  REALTYPE, dimension(:,:,:), allocatable   :: hcc
  REALTYPE, dimension(:,:,:), allocatable   :: uuEx,vvEx
  REALTYPE, dimension(:,:,:), allocatable   :: num,nuh
  REALTYPE, dimension(:,:,:), allocatable   :: tke,eps
  REALTYPE, dimension(:,:,:), allocatable   :: SS

#ifndef NO_BAROCLINIC
  REALTYPE, dimension(:,:,:), allocatable   :: NN

! 3D baroclinic fields
  REALTYPE, dimension(:,:,:), allocatable   :: S,T,rho,buoy
  REALTYPE, dimension(:,:,:), allocatable   :: idpdx,idpdy
  REALTYPE, dimension(:,:,:), allocatable   :: rad,light
#endif

! suspended matter
#ifndef NO_SUSP_MATTER
  REALTYPE, dimension(:,:,:), allocatable   :: spm,spm_ws
  REALTYPE, dimension(:,:), allocatable     :: spm_pool
#endif

#ifdef UV_TVD
  REALTYPE, dimension(:,:,:), allocatable   :: uadv,vadv,wadv
  REALTYPE, dimension(:,:,:), allocatable   :: huadv,hvadv,hoadv,hnadv
#endif

! 2D fields in 3D domain
  REALTYPE, dimension(:,:), allocatable     :: sseo,ssen
  REALTYPE, dimension(:,:), allocatable     :: ssuo,ssun
  REALTYPE, dimension(:,:), allocatable     :: ssvo,ssvn

! 3D friction in 3D domain
  REALTYPE, dimension(:,:), allocatable     :: rru,rrv,taus,taub

! attenuation
  REALTYPE, dimension(:,:), allocatable     :: A,g1,g2
  
