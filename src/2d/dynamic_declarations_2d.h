#ifdef USE_BREAKS
  integer,dimension(:,:),allocatable   :: break_mask
  integer,dimension(:,:),allocatable   :: break_stat
#endif
  REALTYPE,dimension(:,:),allocatable,target :: t_zo,t_z,D,Dvel,DU,DV
  REALTYPE,dimension(:,:),allocatable  :: U,V
  REALTYPE,dimension(:,:),allocatable  :: UEx,VEx
  REALTYPE,dimension(:,:),allocatable  :: fU,fV
  REALTYPE,dimension(:,:),allocatable  :: ru,rv
  REALTYPE,dimension(:,:),allocatable  :: Uint,Vint
  REALTYPE,dimension(:,:),allocatable  :: res_du,res_u
  REALTYPE,dimension(:,:),allocatable  :: res_dv,res_v
  REALTYPE,dimension(:,:),allocatable  :: SlUx,SlVx
  REALTYPE,dimension(:,:),allocatable  :: Slru,Slrv
  REALTYPE,dimension(:,:),allocatable  :: fwf,fwf_int

  REALTYPE,dimension(:),  allocatable:: EWbdy,ENbdy,EEbdy,ESbdy

! Remember to update this value if you add more 2D arrays.
  integer, parameter :: n2d_fields=27

