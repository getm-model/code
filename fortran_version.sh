#!/bin/sh

if [ "x$FORTRAN_COMPILER" = "xIFORT" ] ; then 
   ifort -v 2> ifort.tmp
   echo  "#define FORTRAN_VERSION \"`head ifort.tmp`\""
   rm -f ifort.tmp
fi

if [ "x$FORTRAN_COMPILER" = "xGFORTRAN" ] ; then 
   echo  "#define FORTRAN_VERSION \"gfortran `gfortran -dumpversion`\""
fi
