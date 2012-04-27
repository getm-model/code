#!/bin/sh

if [ "x$FORTRAN_COMPILER" = "xIFORT" ] ; then 
   ifort --version 2>&1 > ifort.tmp
   echo  "#define FORTRAN_VERSION \"`head -1 ifort.tmp`\""
   rm -f ifort.tmp
fi

if [ "x$FORTRAN_COMPILER" = "xGFORTRAN" ] ; then 
  gfortran --version 2>&1 >  gfortran.tmp
  echo  "#define FORTRAN_VERSION \"gfortran `head -1 gfortran.tmp`\""
  rm gfortran.tmp
fi
