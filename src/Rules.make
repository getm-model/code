#$Id: Rules.make,v 1.4 2003-04-01 15:33:19 gotm Exp $
#
# This file contains rules which are shared between multiple Makefiles.
# This file is quite complicated - all compilation options are set in this
# file and included in all other Makefiles.
# Here you set the fortran compiler and the compilation mode (debug,prof,prod)
# Presently, reasonable values for the NAG Fortran 95 and the Compaq Fortran
# compilers are included.
# If you want to include info for a new compiler - please know what you are 
# doing.
#

SHELL   = /bin/sh

VERSION = 1
PATCHLEVEL = 1
SUBLEVEL = 0

VER     = $(VERSION).$(PATCHLEVEL).$(SUBLEVEL)

# The compilation mode is obtained from $COMPILATION_MODE
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

# The compilation mode is obtained from $COMPILATION_MODE
# default production - else debug or profiling
ifndef GETM_PARALLEL
parallel=false
else
parallel=true
endif

turbulence=
turbulence=gotm

CPP	= /lib/cpp

# Here you can put defines for the [c|f]pp - some will also be set depending
# on compilation mode - if STATIC is defined be careful.
DEFINES =
#DEFINES += -DSPHERICAL
#DEFINES += -DCURVILINEAR
#DEFINES += -DNO_BOTTFRIC
#DEFINES += -DNO_ADVECT
#DEFINES += -DNO_SLR
#DEFINES += -DNEW_CORI
#DEFINES += -DCONST_VISC
#DEFINES += -DPARABOLIC_VISCOSITY
#DEFINES += -DNEW_METHOD_KBK
#DEFINES += -DMIN_VEL_DEPTH
#DEFINES += -DPRESS_GRAD_Z
#DEFINES += -DITERATE_VERT_ADV
#DEFINES += -DSUBSTR_INI_PRESS
#DEFINES += -DSONG_WRIGHT
#DEFINES += -DUV_TVD
#DEFINES += -ECMWF_FRV
#DEFINES += -HIRLAM_FRV
#DEFINES += -SETTING_LON_LAN
ifdef STATIC
else
endif

# Directory related settings.

# Top of this version of getm.
ifndef GETMDIR
GETMDIR  = $(HOME)/getm
endif

ifndef BINDIR
BINDIR	= $(GETMDIR)/bin
endif

ifndef LIBDIR
LIBDIR  = $(GETMDIR)/lib/$(FORTRAN_COMPILER)
endif

ifndef MODDIR
MODDIR	= $(GETMDIR)/modules/$(FORTRAN_COMPILER)
endif

INCDIRS		= -I$(GETMDIR)/include -I$(MODDIR)
LINKDIRS	= -L$(LIBDIR)
EXTRA_LIBS	=

ifeq ($(turbulence),gotm)
ifndef GOTMDIR
GOTMDIR	= $(HOME)/gotm
endif
GOTMLIBDIR	= $(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LINKDIRS	+= -L$(GOTMLIBDIR)
EXTRA_LIBS	+= -lturbulence$(buildtype) -lutil$(buildtype) 
INCDIRS		+= -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)
else
DEFINES		+= -DOLD_TURBULENCE
EXTRA_LIBS	+= -lturb$(buildtype)
endif

# Where does the NetCDF include file and library reside.
ifdef NETCDFINC
INCDIRS		+= -I$(NETCDFINC)
endif
ifdef NETCDFLIBNAME
NETCDFLIB	= $(NETCDFLIBNAME)
else
NETCDFLIB	= -lnetcdf
ifdef NETCDFLIBDIR
LDFLAGS		+= -L$(NETCDFLIBDIR)
endif
endif
EXTRA_LIBS	+= $(NETCDFLIB)


# Where does the MPI library reside.
ifeq ($(parallel),true)
DEFINES += -DPARALLEL
ifdef MPIINC
INCDIRS		+= -I$(MPIINC)
endif
ifdef MPIMOD
INCDIRS		+= -I$(MPIMOD)
endif
ifdef MPILIBNAME
MPILIB		= $(MPILIBNAME)
else
MPILIB		= -lmpich
MPILIB		=
ifdef MPILIBDIR
LDFLAGS		+= -L$(MPILIBDIR)
LINKDIRS	+= -L$(MPILIBDIR)
endif
endif
EXTRA_LIBS	+= $(MPILIB)
endif

DOCDIR		= $(GETMDIR)/doc

# Information about which Fortran compiler to use is 
# obtained from $(FORTRAN_COMPILER) - environment variable.

# Normaly this should not be changed - unless you want something very specific.

# The Fortran compiler is determined from the EV FORTRAN_COMPILER - options 
# sofar NAG(linux), FUJITSU(Linux), DECF90 (OSF1 and likely Linux on alpha),
# SunOS, PGF90 - Portland Group Fortran Compiler (on Intel Linux).

# Set options for the NAG Fortran compiler.
ifeq ($(FORTRAN_COMPILER),NAG)
FC=f95nag
DEFINES += -DFORTRAN95
can_do_F90=true
MODULES=-mdir $(MODDIR)
EXTRAS	= -f77
DEBUG_FLAGS = -g -C=all -O0
PROF_FLAGS  = -pg -O3
PROD_FLAGS  = -O3
REAL_4B	= real*4
endif


# Set options for the Compaq fort compiler - on alphas.
ifeq ($(FORTRAN_COMPILER),DECFOR)
FC=f90
DEFINES += -DFORTRAN95
can_do_F90=false
MODULES=-module $(MODDIR)
EXTRAS	=
DEBUG_FLAGS = -g -arch host -check bounds -check overflow -check nopower -check underflow -std90 -assume gfullpath 
DEBUG_FLAGS = -g -arch host -check bounds -check overflow -check nopower -assume gfullpath 
PROF_FLAGS  = -pg -O
PROD_FLAGS  = -O -fast -inline speed -pipeline
#PROD_FLAGS  = -O -fast -inline speed -unroll  1 -pipeline
REAL_4B	= real\(4\)
endif

# Set options for the Fujitsu compiler - on Linux/Intel and SunOS.
ifeq ($(FORTRAN_COMPILER),FUJITSU)
FC=frt
can_do_F90=true
DEFINES += -DFORTRAN95
MODULES=-Am -M$(MODDIR)
EXTRAS  = -ml=cdecl -fw
EXTRAS  = -fw
DEBUG_FLAGS = -g -H aeus
PROF_FLAGS  = -pg -O3
PROD_FLAGS  = -O -K fast
REAL_4B	= real\(4\)
endif

# Set options for the Portland Group Fortran 90 compiler.
ifeq ($(FORTRAN_COMPILER),PGF90)
FC=pgf90
DEFINES += -DFORTRAN90
can_do_F90=false
can_do_F90=true
F90_to_f90=$(FC) -E $(F90FLAGS) $(EXTRA_FFLAGS) $< > $@
MODULES=-module $(MODDIR)
EXTRAS  =
DEBUG_FLAGS = -g
PROF_FLAGS  =
PROD_FLAGS  = -fast
REAL_4B = real\(4\)
endif

# Set options for the Intel Fortran 95 compiler.
ifeq ($(FORTRAN_COMPILER),IFC)
FC=ifc
DEFINES += -DFORTRAN95
can_do_F90=true
F90_to_f90=$(FC) -E $(F90FLAGS) $(EXTRA_FFLAGS) $< > $@
F90_to_f90=
MODULES=
MODULES=-module $(MODDIR)
EXTRAS  = -static -w95 -e95
DEBUG_FLAGS = -g -C
PROF_FLAGS  = -qp -p
PROD_FLAGS  = -O3 -mp
REAL_4B = real\(4\)
EXTRA_LIBS += -lPEPCF90 -lpthread
endif

DEFINES += -DREAL_4B=$(REAL_4B)

# Sets options for debug compilation
ifeq ($(compilation),debug)
buildtype = _debug
DEFINES += -DDEBUG $(STATIC)
FLAGS   = $(DEBUG_FLAGS) 
endif

# Sets options for profiling compilation
ifeq ($(compilation),profiling)
buildtype = _prof
DEFINES += -DPROFILING $(STATIC)
FLAGS   = $(PROF_FLAGS) 
endif

# Sets options for production compilation
ifeq ($(compilation),production)
buildtype = _prod
DEFINES += -DPRODUCTION $(STATIC)
FLAGS   = $(PROD_FLAGS) 
endif

# For making the source code documentation.
PROTEX	= protex -b -n -s

# It should not be necessary to change anything below this line - kbk

#
# False targets.
#
.PHONY: dummy

.SUFFIXES:
.SUFFIXES: .F90


CPPFLAGS	= $(DEFINES) $(INCDIRS)
FFLAGS  	= $(DEFINES) $(FLAGS) $(MODULES) $(INCDIRS) $(EXTRAS)
F90FLAGS  	= $(FFLAGS)
LDFLAGS		= $(FFLAGS) $(LINKDIRS)

#
# Special variables which should not be exported
#
unexport SUBDIRS

#
# Common rules
#

ifeq  ($(can_do_F90),true)
%.o: %.F90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $<
else
%.f90: %.F90
	$(F90_to_f90)
#	$(CPP) $(CPPFLAGS) $< -o $@

%.o: %.f90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
endif

