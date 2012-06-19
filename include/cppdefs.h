! This file is include in all .F90 files and contains very important
! difinitions. Infact the model will not compile when this file is not
! in a correct format.
! KBK 990615

#include "version.h"
#include "git_revision.h"
#include "fortran_version.h"

#if defined(SPHERICAL) || defined(CURVILINEAR)
#define DXC dxc(i,j)
#define DXCJP1 dxc(i,j+1)
#define DXU dxu(i,j)
#define DXUIM1 dxu(i-1,j)
#define DXUJP1 dxu(i,j+1)
#define DXV dxv(i,j)
#define DXVIP1 dxv(i+1,j)
#define DXVJM1 dxv(i,j-1)
#define DXVJP1 dxv(i,j+1)
#define DXVPM dxv(i+1,j-1)
#define DXX dxx(i,j)
#define DXXJM1 dxx(i,j-1)
#define DYC dyc(i,j)
#define DYCIP1 dyc(i+1,j)
#define DYU dyu(i,j)
#define DYUIP1 dyu(i+1,j)
#define DYUIM1 dyu(i-1,j)
#define DYUJP1 dyu(i,j+1)
#define DYUMP dyu(i-1,j+1)
#define DYV dyv(i,j)
#define DYVIP1 dyv(i+1,j)
#define DYVJM1 dyv(i,j-1)
#define DYX dyx(i,j)
#define DYXIM1 dyx(i-1,j)
#define ARCD1 arcd1(i,j)
#define ARUD1 arud1(i,j)
#define ARVD1 arvd1(i,j)
#else
#define DXC dx
#define DXCJP1 dx
#define DXU dx
#define DXUIM1 dx
#define DXUJP1 dx
#define DXV dx
#define DXVIP1 dx
#define DXVJM1 dx
#define DXVJP1 dx
#define DXVPM dx
#define DXX dx
#define DXXJM1 dx
#define DYC dy
#define DYCIP1 dy
#define DYU dy
#define DYUIP1 dy
#define DYUIM1 dy
#define DYUJP1 dy
#define DYUMP dy
#define DYV dy
#define DYVIP1 dy
#define DYVJM1 dy
#define DYX dy
#define DYXIM1 dy
#define ARCD1 ard1
#define ARUD1 ard1
#define ARVD1 ard1
#endif

! For 2D boundary conditions
#define ZERO_GRADIENT 1
#define SOMMERFELD    2
#define CLAMPED_ELEV  3
#define FLATHER_ELEV  4
#define FLATHER_VEL   5
#define CLAMPED_VEL   6
#define CLAMPED       7

! Reserved Fortran units
#define stdin  		5
#define stdout 		6
#define stderr 		0
#define debug  		0
#define NAMLST 		10
#define NAMLST2		11
#define ADAPTNML	12
#define RESTART 	15
#define PARSETUP 	20
#define BDYINFO 	21
#define BDYDATA 	22

! Data/file formats
#define NO_DATA		-1
#define ANALYTICAL	0
#define ASCII		1
#define NETCDF		2
#define BINARY		3
#define OPENDX		4
#define GRADS		5

#define PATH_MAX	255

! Handy for writing
#define STDOUT write(stdout,*)
#define STDERR write(stderr,*)
#define LEVEL0 STDERR
#define LEVEL1 STDERR '   ',
#define LEVEL2 STDERR '       ',
#define LEVEL3 STDERR '           ',
#define LEVEL4 STDERR '               ',
#define FATAL  STDERR 'FATAL ERROR: ',

!KBK#define STDERR IF(myid.le.0) write(stderr,*)
#define LINE "------------------------------------------------------------------------"

! For easier reading
#define READING 0
#define WRITING 1

! To avoid dividing by zero
#define SMALL 1e-8

! What precision will we use in this compilation
#define SINGLE
#undef  SINGLE

#ifdef SINGLE
#define REAL_SIZE 4
#define REALTYPE real(kind=4)
#define _ZERO_  0.0
#define _TENTH_ 0.1
#define _QUART_ 0.25
#define _HALF_  0.5
#define _ONE_   1.0
#define _TWO_   2.0
#define _THREE_ 3.0
#else
#define REAL_SIZE 8
#define REALTYPE real(kind=8)
#define MPI_REALTYPE MPI_DOUBLE_PRECISION
#define _ZERO_  0.0d0
#define _TENTH_ 0.1d0
#define _QUART_ 0.25d0
#define _HALF_  0.5d0
#define _ONE_   1.0d0
#define _TWO_   2.0d0
#define _THREE_ 3.0d0
#endif

! Rare use of long integers (e.g. timers):
#define LONGINT INTEGER(KIND=8)

! The width of the HALO zones
#define HALO	  2

! Size of dimensions with and without HALO-zones:
#define _IRANGE_HALO_ imin-HALO:imax+HALO
#define _JRANGE_HALO_ jmin-HALO:jmax+HALO
#define _IRANGE_NO_HALO_ imin:imax
#define _JRANGE_NO_HALO_ jmin:jmax
#define _KRANGE_ 0:kmax

! Here the memory-allocation is defined using above definitions
#define E2DFIELD  _IRANGE_HALO_,_JRANGE_HALO_
#define E2DXFIELD -1+_IRANGE_HALO_,-1+_JRANGE_HALO_
#define I2DFIELD  _IRANGE_HALO_,_JRANGE_HALO_
#define I3DFIELD  _IRANGE_HALO_,_JRANGE_HALO_,_KRANGE_

! _2D_W_HOT_, _3D_W_HOT_ - macros for writing hot-start files
#ifdef _WRITE_HOT_HALOS_
#warning "_WRITE_HOT_HALOS_ not implemented yet - undef it"
#undef _WRITE_HOT_HALOS_
#endif
#define _2D_W_HOT_ _IRANGE_NO_HALO_,_JRANGE_NO_HALO_
#define _3D_W_HOT_ _IRANGE_NO_HALO_,_JRANGE_NO_HALO_,_KRANGE_

! _2D_R_HOT_, _3D_R_HOT_ - macros for reading hot-start files
#ifdef _READ_HOT_HALOS_
#warning "READ_HOT_HALOS_ not implemented yet - undef it"
#undef _READ_HOT_HALOS_
#endif
#define _2D_R_HOT_ _IRANGE_NO_HALO_,_JRANGE_NO_HALO_
#define _3D_R_HOT_ _IRANGE_NO_HALO_,_JRANGE_NO_HALO_,_KRANGE_

! _2D_W_, _3D_W_ to specify array slices in nf90_put_var()
#ifdef _WRITE_HALOS_
#warning "_WRITE_HALOS_ not implemented yet - undef it"
#undef _WRITE_HALOS_
! need to update/fix .../futils/{cnv_,to}*.F90
! loop boundaries and mask issue
#endif

#ifdef _WRITE_HALOS_
#define _IRANGE_ _IRANGE_HALO_
#define _JRANGE_ _JRANGE_HALO_
#else
#define _IRANGE_ _IRANGE_NO_HALO_
#define _JRANGE_ _JRANGE_NO_HALO_
#endif
#define _2D_W_ _IRANGE_,_JRANGE_
#define _3D_W_ _IRANGE_,_JRANGE_,_KRANGE_

! Definition of vertical coordinate identifiers
#define _SIGMA_COORDS_        1
#define _Z_COORDS_            2
#define _GENERAL_COORDS_      3
#define _HYBRID_COORDS_       4
#define _ADAPTIVE_COORDS_     5

! Definition to write NetCDF output reals as single or double precision:
#ifdef _NCDF_SAVE_DOUBLE_
#define NCDF_FLOAT_PRECISION NF90_DOUBLE
#else
#define NCDF_FLOAT_PRECISION NF90_REAL
#endif
