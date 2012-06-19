! Remember to update this value if you add more 3D arrays.
  integer,parameter                    :: n3d_fields=29

! Number of vertical layers in z,u,v columns
  integer, dimension(:,:), allocatable:: kmin,kumin,kvmin
  integer, dimension(:,:), allocatable:: kmin_pmz,kumin_pmz,kvmin_pmz

  REALTYPE, dimension(:,:,:), allocatable   :: uu,vv
  REALTYPE, dimension(:,:,:), allocatable, target :: ww

#ifdef _MOMENTUM_TERMS_
  REALTYPE, dimension(:,:,:), allocatable, target :: tdv_u,adv_u,vsd_u
  REALTYPE, dimension(:,:,:), allocatable, target :: hsd_u,cor_u,epg_u
  REALTYPE, dimension(:,:,:), allocatable, target :: ipg_u

  REALTYPE, dimension(:,:,:), allocatable, target :: tdv_v,adv_v,vsd_v
  REALTYPE, dimension(:,:,:), allocatable, target :: hsd_v,cor_v,epg_v
  REALTYPE, dimension(:,:,:), allocatable, target :: ipg_v
#endif

#ifdef STRUCTURE_FRICTION
  REALTYPE, dimension(:,:,:), allocatable   :: sf
#endif
  REALTYPE, dimension(:,:,:), allocatable, target :: hn,hun,hvn
  REALTYPE, dimension(:,:,:), allocatable   :: ho,huo,hvo
  REALTYPE, dimension(:,:,:), allocatable   :: hcc
  REALTYPE, dimension(:,:,:), allocatable   :: uuEx,vvEx
  REALTYPE, dimension(:,:,:), allocatable, target :: nuh
  REALTYPE, dimension(:,:,:), allocatable   :: num
  REALTYPE, dimension(:,:,:), allocatable   :: tke,eps
  REALTYPE, dimension(:,:,:), allocatable   :: SS

#ifndef NO_BAROCLINIC
  REALTYPE, dimension(:,:,:), allocatable   :: NN

! 3D baroclinic fields
  REALTYPE, dimension(:,:,:), allocatable, target :: S,T,rho
  REALTYPE, dimension(:,:,:), allocatable   :: buoy
  REALTYPE, dimension(:,:,:), allocatable   :: alpha,beta
  REALTYPE, dimension(:,:,:), allocatable   :: idpdx,idpdy
  REALTYPE, dimension(:,:,:), allocatable   :: rad,light
  REALTYPE, dimension(:,:,:), allocatable   :: nummix3d_S,nummix3d_T
  REALTYPE, dimension(:,:,:), allocatable   :: phymix3d_S,phymix3d_T
  REALTYPE, dimension(:,:), allocatable     :: nummix2d_S,nummix2d_T
  REALTYPE, dimension(:,:), allocatable     :: phymix2d_S,phymix2d_T
#endif
  REALTYPE, dimension(:,:,:), allocatable   :: numdis3d
  REALTYPE, dimension(:,:), allocatable     :: numdis2d

! suspended matter
#ifndef NO_SUSP_MATTER
  REALTYPE, dimension(:,:,:), allocatable   :: spm,spm_ws
  REALTYPE, dimension(:,:), allocatable     :: spm_pool
#endif

! 2D fields in 3D domain
  REALTYPE, dimension(:,:), allocatable     :: sseo,ssen,Dn
  REALTYPE, dimension(:,:), allocatable     :: ssuo,ssun
  REALTYPE, dimension(:,:), allocatable     :: ssvo,ssvn
  REALTYPE,dimension(:,:),allocatable,target :: Dun,Dvn

! 3D friction in 3D domain
  REALTYPE, dimension(:,:), allocatable     :: rru,rrv,zub,zvb
  REALTYPE, dimension(:,:), allocatable     :: taus,taubx,tauby,taub

! attenuation
  REALTYPE, dimension(:,:), allocatable     :: A,g1,g2

