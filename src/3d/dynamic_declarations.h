! Remember to update this value if you add more 3D arrays.
#ifdef UV_TVD
  integer, parameter ::	n3d_fields=29
#else
  integer, parameter ::	n3d_fields=22
#endif
  REALTYPE, dimension(:,:,:), allocatable  	:: uu,vv,ww
  REALTYPE, dimension(:,:,:), allocatable	:: ho,hn
  REALTYPE, dimension(:,:,:), allocatable	:: huo,hun
  REALTYPE, dimension(:,:,:), allocatable	:: hvo,hvn
  REALTYPE, dimension(:,:,:), allocatable	:: uuEx,vvEx
  REALTYPE, dimension(:,:,:), allocatable	:: num,nuh
  REALTYPE, dimension(:,:,:), allocatable	:: tke,eps
#ifdef OLD_TURBULENCE
  REALTYPE, dimension(:,:,:), allocatable	:: tkeo,P,B
#endif
  REALTYPE, dimension(:,:,:), allocatable	:: SS,NN

! 3D baroclinic fields
  REALTYPE, dimension(:,:,:), allocatable	:: S,T,rho
  REALTYPE, dimension(:,:,:), allocatable	:: idpdx,idpdy

! suspended matter
  REALTYPE, dimension(:,:,:), allocatable	:: spm,spm_ws
  REALTYPE, dimension(:,:), allocatable		:: spm_pool

#ifdef UV_TVD
  REALTYPE, dimension(:,:,:), allocatable	:: uadv,vadv,wadv
  REALTYPE, dimension(:,:,:), allocatable	:: huadv,hvadv,hoadv,hnadv
#endif

! 2D fields in 3D domain
  REALTYPE, dimension(:,:), allocatable		:: sseo,ssen
  REALTYPE, dimension(:,:), allocatable		:: ssuo,ssun
  REALTYPE, dimension(:,:), allocatable		:: ssvo,ssvn

! 3D friction in 3D domain
  REALTYPE, dimension(:,:), allocatable		:: rru,rrv,taus,taub

! Number of vertical layers in z,u,v columns
  integer, dimension(:,:), allocatable		:: kmin,kumin,kvmin
