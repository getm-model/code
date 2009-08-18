#$Id: Rules.make,v 1.21 2009-08-18 10:24:43 bjb Exp $
#
# This file contains rules which are shared between multiple Makefiles.
# This file is quite complicated - all compilation options are set in this
# file and included in all other Makefiles.
# Here you set the fortran compiler and the compilation mode (debug,prof,prod)
# If you want to include info for a new compiler - please know what you are 
# doing - then add a new file to the directory ../compilers.

SHELL   = /bin/sh
CPP	= /lib/cpp

# Information about which Fortran compiler to use is 
# obtained from $(FORTRAN_COMPILER) - environment variable.
# The file ../compilers/compiler.$(FORTRAN_COMPILER) must exist

DEFINES=-D$(FORTRAN_COMPILER)
include $(GETMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

# The compilation mode is obtained from $COMPILATION_MODE
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

# 2D compilation only
ifeq ($(GETM_NO_3D),true)
DEFINES += -DNO_3D
export GETM_NO_BAROCLINIC=true
endif

# 3D barotropic
ifeq ($(GETM_NO_BAROCLINIC),true)
DEFINES += -DNO_BAROCLINIC
endif

# Suspended matter
ifeq ($(GETM_SPM),true)
DEFINES += -DSPM
endif

# Structure friction
ifeq ($(GETM_STRUCTURE_FRICTION),true)
DEFINES += -DSTRUCTURE_FRICTION
endif

# Bio-geochemical component
ifeq ($(GETM_BIO),true)
DEFINES += -DGETM_BIO
endif

# Remove timers
ifeq ($(GETM_NO_TIMERS),true)
DEFINES += -DNO_TIMERS
endif

# Compile for parallel execution
ifeq ($(GETM_PARALLEL),true)
parallel=true
set par=par
else
parallel=false
set par=ser
endif

turbulence=
turbulence=gotm

# Here you can put defines for the [c|f]pp - some will also be set depending
# on compilation mode - if STATIC is defined be careful.
ifdef INPUT_DIR
DEFINES += -DINPUT_DIR="'$(INPUT_DIR)/'"
endif
#DEFINES += -DSPHERICAL
#DEFINES += -DCURVILINEAR
#DEFINES += -DNO_BOTTFRIC
#DEFINES += -DNO_ADVECT
#DEFINES += -DNO_SLR
#DEFINES += -DNEW_CORI
#DEFINES += -DCONSTANT_VISCOSITY
#DEFINES += -DPARABOLIC_VISCOSITY
#DEFINES += -DNEW_METHOD_KBK
#DEFINES += -DMIN_VEL_DEPTH
#DEFINES += -DPRESS_GRAD_Z
#DEFINES += -DITERATE_VERT_ADV
#DEFINES += -DSUBSTR_INI_PRESS
#DEFINES += -DSONG_WRIGHT
#DEFINES += -DUV_TVD
ifdef STATIC
else
endif

# Directory related settings.

# Top of this version of getm.
ifndef GETMDIR
GETMDIR  = $(HOME)/GETM/getm-devel
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

ifeq ($(GETM_BIO),true)
EXTRA_LIBS += -lbio$(buildtype)
endif

ifeq ($(turbulence),gotm)
ifndef GOTMDIR
GOTMDIR = $(HOME)/gotm
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

ifdef NETCDFLIBDIR
LINKDIRS	+= -L$(NETCDFLIBDIR)
endif

ifdef NETCDFLIBNAME
NETCDFLIB	= $(NETCDFLIBNAME)
else
NETCDFLIB	= -lnetcdf
endif

ifeq ($(NETCDF_VERSION),NETCDF4)

DEFINES		+= -DNETCDF4
ifdef HDF5_DIR
INCDIRS		+= -I$(HDF5_DIR)/include
LINKDIRS	+= -L$(HDF5_DIR)/lib
endif
HDF5LIB		= -lhdf5_hl -lhdf5 -lz

else  # NetCDF3 is default

DEFINES		+= -DNETCDF3
HDF5LIB		=

endif

EXTRA_LIBS	+= $(NETCDFLIB) $(HDF5LIB)

# NetCDF/HDF configuration done

# Where does the MPI library reside.
ifeq ($(parallel),true)
DEFINES += -DPARALLEL

# OPENMPI - set FC to mpif90
ifeq ($(MPI),OPENMPI)
FC=mpif90
endif

# MPICH2 - set FC to mpif90
ifeq ($(MPI),MPICH2)
FC=mpif90
endif

# SGI - MPI - works for Peter Holtermann on ALTIX 3700.
ifeq ($(MPI),SGIMPI)
EXTRA_LIBS      += -lmpi
endif

# obsolete - use either OPenMPI or MPICH2
ifeq ($(MPI),MPICH)
ifdef MPIINC
INCDIRS		+= -I$(MPIINC)
endif
ifdef MPIMOD
INCDIRS		+= -I$(MPIMOD)
endif
ifdef MPILIBFILE
MPILIB		= $(MPILIBFILE)
else
ifdef MPILIBNAME
MPILIB		= $(MPILIBNAME)
else
MPILIB		= -lmpich -lpthread
endif
ifdef MPILIBDIR
LDFLAGS		+= -L$(MPILIBDIR)
LINKDIRS	+= -L$(MPILIBDIR)
endif
endif
EXTRA_LIBS	+= $(MPILIB)
endif

endif

DOCDIR		= $(GETMDIR)/doc

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
