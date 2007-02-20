! Remember to update this value if you add more 3D arrays.
#ifdef UV_TVD
#ifdef SPM
  integer, parameter                   :: n3d_fields=29
#else
  integer, parameter                   :: n3d_fields=26 
#endif  
#else
#ifdef SPM
  integer, parameter                   :: n3d_fields=22
#else
  integer, parameter                   :: n3d_fields=19
#endif
#endif
! Number of vertical layers in z,u,v columns
  INTEGER                              :: kmin(I2DFIELD)
  INTEGER                              :: kumin(I2DFIELD)
  INTEGER                              :: kvmin(I2DFIELD)
  INTEGER                              :: kmin_pmz(I2DFIELD)
  INTEGER                              :: kumin_pmz(I2DFIELD)
  INTEGER                              :: kvmin_pmz(I2DFIELD)

  REALTYPE                             :: uu(I3DFIELD)
  REALTYPE                             :: vv(I3DFIELD)
  REALTYPE                             :: ww(I3DFIELD)
  REALTYPE                             :: ho(I3DFIELD)
  REALTYPE                             :: hn(I3DFIELD)
  REALTYPE                             :: huo(I3DFIELD)
  REALTYPE                             :: hun(I3DFIELD)
  REALTYPE                             :: hvo(I3DFIELD)
  REALTYPE                             :: hvn(I3DFIELD)
  REALTYPE                             :: hcc(I3DFIELD)
  REALTYPE                             :: uuEx(I3DFIELD)
  REALTYPE                             :: vvEx(I3DFIELD)
  REALTYPE                             :: num(I3DFIELD)
  REALTYPE                             :: nuh(I3DFIELD)

! 3D turbulent fields
  REALTYPE                             :: tke(I3DFIELD)
  REALTYPE                             :: eps(I3DFIELD)
  REALTYPE                             :: SS(I3DFIELD)

#ifndef NO_BAROCLINIC
! 3D baroclinic fields
  REALTYPE                             :: NN(I3DFIELD)
  REALTYPE                             :: S(I3DFIELD)
  REALTYPE                             :: T(I3DFIELD)
  REALTYPE                             :: rho(I3DFIELD)
  REALTYPE                             :: idpdx(I3DFIELD)
  REALTYPE                             :: idpdy(I3DFIELD)
  REALTYPE                             :: rad(I3DFIELD)
#endif

#ifdef SPM
! suspended matter
  REALTYPE                             :: spm(I3DFIELD)
  REALTYPE                             :: spm_ws(I3DFIELD)
  REALTYPE                             :: spm_pool(I2DFIELD)
#endif

  REALTYPE                             :: light(I3DFIELD)

#ifdef UV_TVD
  REALTYPE                             :: uadv(I3DFIELD)
  REALTYPE                             :: vadv(I3DFIELD)
  REALTYPE                             :: wadv(I3DFIELD)
  REALTYPE                             :: huadv(I3DFIELD)
  REALTYPE                             :: hvadv(I3DFIELD)
  REALTYPE                             :: hoadv(I3DFIELD)
  REALTYPE                             :: hnadv(I3DFIELD)
#endif

! 2D fields in 3D domain
  REALTYPE                             :: sseo(I2DFIELD)
  REALTYPE                             :: ssen(I2DFIELD)
  REALTYPE                             :: ssuo(I2DFIELD)
  REALTYPE                             :: ssun(I2DFIELD)
  REALTYPE                             :: ssvo(I2DFIELD)
  REALTYPE                             :: ssvn(I2DFIELD)

! 3D friction in 3D domain
  REALTYPE                             :: rru(I2DFIELD)
  REALTYPE                             :: rrv(I2DFIELD)
  REALTYPE                             :: taus(I2DFIELD)
  REALTYPE                             :: taub(I2DFIELD)

! light attenuation
  REALTYPE                             :: A(I2DFIELD)
  REALTYPE                             :: g1(I2DFIELD)
  REALTYPE                             :: g2(I2DFIELD)

