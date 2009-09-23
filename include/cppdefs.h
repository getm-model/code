! This file is include in all .F90 files and contains very important
! difinitions. Infact the model will not compile when this file is not
! in a correct format.
! KBK 990615

#include "version.h"

#if defined(SPHERICAL) || defined(CURVILINEAR)
#define DXC dxc(i,j)
#define DXCIP1 dxc(i+1,j)
#define DXCIM1 dxc(i-1,j)
#define DXCJP1 dxc(i,j+1)
#define DXU dxu(i,j)
#define DXV dxv(i,j)
#define DXVJM1 dxv(i,j-1)
#define DXX dxx(i,j)
#define DXXJM1 dxx(i,j-1)
#define DYC dyc(i,j)
#define DYCIP1 dyc(i+1,j)
#define DYCJP1 dyc(i,j+1)
#define DYCJM1 dyc(i,j-1)
#define DYU dyu(i,j)
#define DYUIM1 dyu(i-1,j)
#define DYV dyv(i,j)
#define DYX dyx(i,j)
#define DYXIM1 dyx(i-1,j)
#define ARCD1 arcd1(i,j)
#define ARUD1 arud1(i,j)
#define ARVD1 arvd1(i,j)
#else
#define DXC dx
#define DXCIP1 dx
#define DXCIM1 dx
#define DXCJP1 dx
#define DXU dx
#define DXV dx
#define DXVJM1 dx
#define DXX dx
#define DXXJM1 dx
#define DYC dy
#define DYCIP1 dy
#define DYCJP1 dy
#define DYCJM1 dy
#define DYU dy
#define DYUIM1 dy
#define DYV dy
#define DYX dy
#define DYXIM1 dy
#define ARCD1 ard1
#define ARUD1 ard1
#define ARVD1 ard1
#endif

! For 2D boundary conditions
#define ZERO_GRADIENT 1
#define SOMMERFELDT   2
#define CLAMPED       3
#define FLATHER_ELEV  4
#define FLATHER_VEL   5

! Reserved Fortran units
#define stdin  		5
#define stdout 		6
#define stderr 		0
#define debug  		0
#define NAMLST 		10
#define NAMLST2		11
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
#define REALTYPE REAL
#define REAL_SIZE 4
#define _ZERO_ 0.0
#define _TENTH_ 0.1
#define _QUART_ 0.25
#define _HALF_ 0.5
#define _ONE_  1.0
#define _TWO_  2.0
#else
#define REALTYPE DOUBLE PRECISION
#define MPI_REALTYPE MPI_DOUBLE_PRECISION
#define REAL_SIZE 8
#define _ZERO_ 0.0d0
#define _TENTH_ 0.1d0
#define _QUART_ 0.25d0
#define _HALF_ 0.5d0
#define _ONE_  1.0d0
#define _TWO_  2.0d0
#endif

! Rare use of long integers (e.g. timers):
#define LONGINT INTEGER(KIND=8)

! The width of the HALO zones
#define HALO	  2

! Here the memory-allocation is defined
#define E2DFIELD  imin-HALO:imax+HALO,jmin-HALO:jmax+HALO
#define E2DXFIELD imin-1-HALO:imax+HALO,jmin-1-HALO:jmax+HALO
#define I2DFIELD  imin-HALO:imax+HALO,jmin-HALO:jmax+HALO
#define I3DFIELD  imin-HALO:imax+HALO,jmin-HALO:jmax+HALO,0:kmax

! For easier saving HALO zones in NetCDF files - used in:
! save_grid_ncdf.F90 and write_restart_ncdf.F90
#ifdef SAVE_HALOS
#define IRANGE imin-HALO:imax+HALO
#define JRANGE jmin-HALO:jmax+HALO
#else
#define IRANGE imin:imax
#define JRANGE jmin:jmax
#endif


! These defines the do loops for the real inner points..
! that is the points that are independent of neighbours.
#define DO_EILOOP  DO i=imin,imax
#define DO_EJLOOP  DO j=jmin,jmax
