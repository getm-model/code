#!/bin/sh

# if not set use the suggested source code installation directories
GETM_BASE=${GETM_BASE:=~/source/repos/GETM/code}

# default Fortran compiler is gfortran - overide by setting compiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}

# horizontal coordinate system to use - default Cartesian
# other options are Spherical or Curvilinear
# to set use e.g.:
# export coordinates=Spherical
coordinate=${coordinate:=Cartesian}

# configurable installation prefix
# override by e.g.:
# export install_prefix=/tmp
# note that $compiler will be appended
install_prefix=${install_prefix:=~/local/getm}

# ready to configure
mkdir -p $compiler
cd $compiler
cmake $GETM_BASE \
      -DGETM_EMBED_VERSION=on \
      -DGETM_USE_FABM=on \
      -DCMAKE_Fortran_COMPILER=$compiler \
      -DGETM_USE_PARALLEL=off \
      -DGETM_COORDINATE_TYPE=$coordinate \
      -DCMAKE_INSTALL_PREFIX=$install_prefix/$compiler
cd ..
