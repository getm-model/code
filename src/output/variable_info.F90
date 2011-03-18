!#ifdef NETCDF4
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !MODULE:  variable_info - metadata for all variables that can be saved.
!
! !INTERFACE:
   module variable_info
!
! !DESCRIPTION:
!  This module will be used by the GETM output system to obtain information
!  about variables. Multiple output definitions i.e. hyperslabs and 
!  variable-lists shall all use the information contained here. The key
!  to the meta-data is the $name$ i.e. search the array until agreement
!  between $varinfo$->$name$ and variable name searched for.
!  This module provides meta data information about all variables GETM 
!  that be saved in output files. The specific data for a given variable
!  is contained in a defined type $varinfo$. $varinfo$ contains a number
!  of fields including a unique $id$, a $name$, information about the 
!  dimensions of the variable and additional fields that can be used by
!  NetCDF to give COARDS conforming meta data information.
!  The $varinfo$ also contains pointer types that can be used to point to
!  the Fortran variable actually holding the data. For a re-useable type
!  definition it is necessary to define pointers to 1D, 2D and 3D 
!  variables.
!  The $varinfo\_list$ can be expanded to also include compound variables
!  i.e. variables that are functions of other variables. In this case the
!  unique $id$ shall be used to actually select which variables and what
!  operations should be done - like:
!  \begin{verbatim}
!  select case (varinfo_list(1)%id)
!     case(?)
!     case(_SALT_FLUX_)
!        data = uu*S
!     case(?)
!  \end{verbatim}
!  \newline
!  The final implementation details not decided yet.\newline
!  Work in progress.
!
! !USES:
   use exceptions
#if 0
!   use domain
!   use variables_2d
!   use variables_3d
!   use meteo
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS
   public init_var_info, print_var_info

   integer, parameter           :: max_length=255
   type, public                 :: varinfo
      integer                   :: id=-1
      character(len=max_length) :: name=''
      integer                   :: ndims=-1
      character(len=max_length) :: long_name=''
      character(len=max_length) :: units=''
      REALTYPE                  :: mv
      REALTYPE                  :: vr(2)
      REALTYPE, pointer         :: data_0d        => null()
      REALTYPE, pointer         :: data_1d(:)     => null()
      REALTYPE, pointer         :: data_2d(:,:)   => null()
      REALTYPE, pointer         :: data_3d(:,:,:) => null()
!      REALTYPE, pointer         :: save_data_2d(:,:)   => null()
   end type varinfo

   type(varinfo)                :: varinfo_list(40)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hannes Rennau
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Optain info on the variable given by var_name
!
! !INTERFACE:
   subroutine init_var_info()
!
! !DESCRIPTION:
!  Here a list meta-data for all possible GETM outputr variables is 
!  constructed. The information is put in an array.
!
   IMPLICIT NONE
#if 1
#define STATIC
#ifdef STATIC
   REALTYPE, target :: z(0:11,-10:10)
#else
#endif
#undef STATIC
#endif

!
! !LOCAL VARIABLES
#define _D_      1
#define _DU_     2
#define _DV_     3
#define _AIRP_   4
#define _U10_    5
#define _V10_    6
#define _T2_     7
#define _HUM_    8
#define _TCC_    9
#define _TAUSX_ 10
#define _TAUSY_ 11
#define _SWR_   12
#define _SHF_   13
#define _ELEV_  14
#define _U_     15
#define _V_     16
#define _RES_U_ 17
#define _RES_V_ 18
#define _UU_    19
#define _VV_    20
#define _HCC_   21
#define _HO_    22
#define _HN_    23
#define _HU_    24
#define _HV_    25
#define _WW_    26
#define _T_     27
#define _S_     28
#define _TKE_   29
#define _NUM_   30
#define _NUH_   31
#define _EPS_   32
#define _SS_    33
#define _NN_    34
#define _SPM_   35
   integer         :: n=0
!EOP
!-----------------------------------------------------------------------
!BOC

!   call random_number(z)
!   n=0

   n=n+1
!  Domain related variables variables
   varinfo_list(n)%id=          _D_
   varinfo_list(n)%name=        'D'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'total water depth'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -15.0, 5000.0 ]
!   varinfo_list(n)%data_2d=>    D

   n=n+1
   varinfo_list(n)%id=          _DU_
   varinfo_list(n)%name=        'DU'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'total water depth - U-points'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -15.0, 5000.0 ]
!   varinfo_list(n)%data_2d=>    DU

   n=n+1
   varinfo_list(n)%id=          _DV_
   varinfo_list(n)%name=        'DV'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'total water depth - V-points'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -15.0, 5000.0 ]
!   varinfo_list(n)%data_2d=>    DV

!  meteo variables
   n=n+1
   varinfo_list(n)%id=          _AIRP_
   varinfo_list(n)%name=        'airp'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'air pressure'
   varinfo_list(n)%units='Pascal'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 90.e3, 110.e3 ]
!   varinfo_list(n)%data_2d=>    airp

   n=n+1
   varinfo_list(n)%id=          _U10_
   varinfo_list(n)%name=        'u10'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'eastward windspeed (10m)'
   varinfo_list(n)%units='m/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -50.0, 50.0 ]
!   varinfo_list(n)%data_2d=>    u10

   n=n+1
   varinfo_list(n)%id=          _V10_
   varinfo_list(n)%name=        'v10'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'northward windspeed (10m)'
   varinfo_list(n)%units='m/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -50.0, 50.0 ]
!   varinfo_list(n)%data_2d=>    v10

   n=n+1
   varinfo_list(n)%id=          _T2_
   varinfo_list(n)%name=        't2'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'temperature (2m)'
   varinfo_list(n)%units='Kelvin'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 325.0 ]
!   varinfo_list(n)%data_2d=>    t2

   n=n+1
   varinfo_list(n)%id=          _HUM_
   varinfo_list(n)%name=        'hum'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'humidity'
   varinfo_list(n)%units='kg/kg'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 100.0 ]
!   varinfo_list(n)%data_2d=>    hum

   n=n+1
   varinfo_list(n)%id=          _TCC_
   varinfo_list(n)%name=        'tcc'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'total cloud cover'
   varinfo_list(n)%units=''
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1.0 ]
!   varinfo_list(n)%data_2d=>    tcc

   n=n+1
   varinfo_list(n)%id=          _TAUSX_
   varinfo_list(n)%name=        'tausx'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'surface stress - x'
   varinfo_list(n)%units='N/m2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -1.0, 1.0 ]
!   varinfo_list(n)%data_2d=>    tausx

   n=n+1
   varinfo_list(n)%id=          _TAUSY_
   varinfo_list(n)%name=        'tausy'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'surface stress - y'
   varinfo_list(n)%units='N/m2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -1.0, 1.0 ]
!   varinfo_list(n)%data_2d=>    tausy

   n=n+1
   varinfo_list(n)%id=          _SWR_
   varinfo_list(n)%name=        'swr'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'short wave radiation'
   varinfo_list(n)%units='W/m2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1500.0 ]
!   varinfo_list(n)%data_2d=>    swr

   n=n+1
   varinfo_list(n)%id=          _SHF_
   varinfo_list(n)%name=        'shf'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'surface heat fluxes'
   varinfo_list(n)%units='W/m2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -1000.0, 1000.0 ]
!   varinfo_list(n)%data_2d=>    shf

!  2D variables
   n=n+1
   varinfo_list(n)%id=          _ELEV_
   varinfo_list(n)%name=        'elev'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'sea surface elevation'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -15.0, 15.0 ]
!   varinfo_list(n)%data_2d=>    elev

   n=n+1
   varinfo_list(n)%id=          _U_
   varinfo_list(n)%name=        'u'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'int. zonal vel.'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3, 3 ]
!   varinfo_list(n)%data_2d=>    u

   n=n+1
   varinfo_list(n)%id=          _V_
   varinfo_list(n)%name=        'v'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'int. meridional vel.'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3, 3 ]
!   varinfo_list(n)%data_2d=>    v

   n=n+1
   varinfo_list(n)%id=          _RES_U_
   varinfo_list(n)%name=        'res_u'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'residual of U'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3, 3 ]
!   varinfo_list(n)%data_2d=>    res_u

   n=n+1
   varinfo_list(n)%id=          _RES_V_
   varinfo_list(n)%name=        'res_v'
   varinfo_list(n)%ndims=       2
   varinfo_list(n)%long_name=   'residual of V'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3, 3 ]
!   varinfo_list(n)%data_2d=>    res_v


!  3D variables
   n=n+1
   varinfo_list(n)%id=          _UU_
   varinfo_list(n)%name=        'uu'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'zonal velocity'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -5.0, 5.0 ]
!   varinfo_list(n)%data_3d=>    uu

   n=n+1
   varinfo_list(n)%id=          _VV_
   varinfo_list(n)%name=        'vv'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'meridional velocity'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -5.0, 5.0 ]
!   varinfo_list(n)%data_3d=>    vv

   n=n+1
   varinfo_list(n)%id=          _HCC_
   varinfo_list(n)%name=        'hcc'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'hydrostatic consistency criterion'
   varinfo_list(n)%units=       ''
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1.0 ]
!   varinfo_list(n)%data_3d=>    hcc

   n=n+1
   varinfo_list(n)%id=          _HO_
   varinfo_list(n)%name=        'ho'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'old box height - h-column'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1000.0 ]
!   varinfo_list(n)%data_3d=>    ho

   n=n+1
   varinfo_list(n)%id=          _HN_
   varinfo_list(n)%name=        'hn'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'new box height - h-column'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1000.0 ]
!   varinfo_list(n)%data_3d=>    hn

   n=n+1
   varinfo_list(n)%id=          _HU_
   varinfo_list(n)%name=        'ho'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'box height - u-column'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1000.0 ]
!   varinfo_list(n)%data_3d=>    hu

   n=n+1
   varinfo_list(n)%id=          _HV_
   varinfo_list(n)%name=        'hv'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'new box height - v-column'
   varinfo_list(n)%units=       'm'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 1000.0 ]
!   varinfo_list(n)%data_3d=>    hv

   n=n+1
   varinfo_list(n)%id=          _WW_
   varinfo_list(n)%name=        'ww'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'vertical velocity'
   varinfo_list(n)%units=       'm/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3.0, 3.0 ]
!   varinfo_list(n)%data_3d=>    ww

   n=n+1
   varinfo_list(n)%id=          _T_
   varinfo_list(n)%name=        'T'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'temperature'
   varinfo_list(n)%units=       'degC'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -2.0, 40.0 ]
!   varinfo_list(n)%data_3d=>    T

   n=n+1
   varinfo_list(n)%id=          _S_
   varinfo_list(n)%name=        'S'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'salinity'
   varinfo_list(n)%units=       'psu'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 40.0 ]
!   varinfo_list(n)%data_3d=>    S

   n=n+1
   varinfo_list(n)%id=          _TKE_
   varinfo_list(n)%name=        'tke'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'turbulent kinetic energy'
   varinfo_list(n)%units=       'm2/s2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 0.2 ]
!   varinfo_list(n)%data_3d=>    tke

   n=n+1
   varinfo_list(n)%id=          _NUM_
   varinfo_list(n)%name=        'num'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'viscosity'
   varinfo_list(n)%units=       'm2/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 0.2 ]
!   varinfo_list(n)%data_3d=>    num

   n=n+1
   varinfo_list(n)%id=          _NUH_
   varinfo_list(n)%name=        'nuh'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'diffusivity'
   varinfo_list(n)%units=       'm2/s'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 0.2 ]
!   varinfo_list(n)%data_3d=>    nuh

   n=n+1
   varinfo_list(n)%id=          _EPS_
   varinfo_list(n)%name=        'eps'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'dissipation'
   varinfo_list(n)%units=       'm2/s3'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ 0.0, 0.2 ]
!   varinfo_list(n)%data_3d=>    eps

   n=n+1
   varinfo_list(n)%id=          _SS_
   varinfo_list(n)%name=        'SS'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'shear-frequency squared'
   varinfo_list(n)%units=       's-2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3.0, 3.0 ]
!   varinfo_list(n)%data_3d=>    SS

   n=n+1
   varinfo_list(n)%id=          _NN_
   varinfo_list(n)%name=        'NN'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'BVF squared'
   varinfo_list(n)%units=       's-2'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3.0, 3.0 ]
!   varinfo_list(n)%data_3d=>    NN

   n=n+1
   varinfo_list(n)%id=          _SPM_
   varinfo_list(n)%name=        'spm'
   varinfo_list(n)%ndims=       3
   varinfo_list(n)%long_name=   'suspended matter concentration'
   varinfo_list(n)%units=       'kg/m-3'
   varinfo_list(n)%mv=          -9999.0
   varinfo_list(n)%vr=          [ -3.0, 3.0 ]
!   varinfo_list(n)%data_3d=>    spm

   return
   end subroutine init_var_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Optainn info on the variable given by var_name
!
! !INTERFACE:
   subroutine print_var_info()
!
! !DESCRIPTION:
!  Just prints out the meta-data for all variables in $varinfo_list$.
!  For testing and debugging purposes mainly.
!
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer          :: n
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'print_var_info() '
   do n=1,size(varinfo_list)
      if (varinfo_list(n)%id .gt. 0) then
      LEVEL3 '----------------------------------------------------------'
      LEVEL3 'variable id:    ',varinfo_list(n)%id
      LEVEL3 'variable name:  ',trim(varinfo_list(n)%name)
      LEVEL3 'number of dims  ',varinfo_list(n)%ndims
      LEVEL3 'long_name:      ',trim(varinfo_list(n)%long_name)
      LEVEL3 'units:          ',trim(varinfo_list(n)%units)
!   LEVEL3 'fill_value=     ',fv
      LEVEL3 'missing_values= ',varinfo_list(n)%mv
      LEVEL3 'lower limit=    ',varinfo_list(n)%vr(1)
      LEVEL3 'upper limit=    ',varinfo_list(n)%vr(2)
      LEVEL3 'data=         ',varinfo_list(n)%data_2d
!      LEVEL3 'save_data=    ',varinfo_list(n)%save_data_2d
      end if
   end do

   do n=1,size(varinfo_list)
      if (varinfo_list(n)%id .gt. 0) then
         select case (varinfo_list(1)%id)
            case(_D_)
               LEVEL3 'var ',n,' is ',trim(varinfo_list(n)%long_name)
            case(_AIRP_)
               LEVEL3 'var ',n,' is ',trim(varinfo_list(n)%long_name)
            case(_ELEV_)
            case(_UU_)
            case default
         end select
      end if
   end do

   return
   end subroutine print_var_info
!EOC

   end module variable_info
!#endif
!-----------------------------------------------------------------------
! Copyright (C) 2009 - Karsten Bolding & Haness Rennau (BB)            !
!-----------------------------------------------------------------------

