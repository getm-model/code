  REALTYPE,dimension(:,:),allocatable  :: D,DU,DV
  REALTYPE,dimension(:,:),allocatable  :: z,zo
  REALTYPE,dimension(:,:),allocatable  :: zu,zv
  REALTYPE,dimension(:,:),allocatable  :: U,V
  REALTYPE,dimension(:,:),allocatable  :: UEx,VEx
  REALTYPE,dimension(:,:),allocatable  :: fU,fV
  REALTYPE,dimension(:,:),allocatable  :: ru,rv
  REALTYPE,dimension(:,:),allocatable  :: Uint,Vint
  REALTYPE,dimension(:,:),allocatable  :: Uinto,Vinto
  REALTYPE,dimension(:,:),allocatable  :: res_du,res_u
  REALTYPE,dimension(:,:),allocatable  :: res_dv,res_v
!kbk
  REALTYPE,dimension(:,:),allocatable  :: ruu,rvv
  REALTYPE,dimension(:,:),allocatable  :: PP
!kbk
  REALTYPE,dimension(:,:),allocatable  :: SlUx,SlVx
  REALTYPE,dimension(:,:),allocatable  :: Slru,Slrv
  REALTYPE,dimension(:,:),allocatable  :: zub,zvb
  REALTYPE,dimension(:,:),allocatable  :: zub0,zvb0
  REALTYPE,dimension(:,:),allocatable  :: surfdiv
  REALTYPE,dimension(:,:),allocatable  :: fwf,fwf_int

  REALTYPE,dimension(:),  allocatable:: EWbdy,ENbdy,EEbdy,ESbdy   

! Remember to update this value if you add more 2D arrays.
  integer, parameter :: n2d_fields=37
