#!/bin/sh

# if not set use the suggested source code installation directories
GETM_BASE=${GETM_BASE:=~/GETM/code}
GOTM_BASE=${GOTM_BASE:=~/GOTM/code}
FABM_BASE=${FABM_BASE:=~/FABM/code}

# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}

# horizontal coordinate system to use - default Cartesian
# other options are Spherical or Curvilinear
# to set use e.g.:
# export coordinates=Spherical
coordinate=${coordinate:=Cartesian}

# NetCDF
# nf-config must be in the path and correpsond to the value of compiler
# try:
# nf-config --all

# ready to configure
mkdir -p $compiler
cd $compiler
cmake $GETM_BASE/src \
      -DGOTM_BASE=$GOTM_BASE \
      -DGETM_USE_FABM=on \
      -DFABM_BASE=$FABM_BASE/ \
      -DCMAKE_Fortran_COMPILER=$compiler \
      -DGETM_USE_PARALLEL=off \
      -DGETM_COORDINATE_TYPE=$coordinate \
      -DCMAKE_INSTALL_PREFIX=~/local/getm/$compiler
cd ..
