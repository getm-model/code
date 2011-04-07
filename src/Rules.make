#$Id: Rules.make,v 1.24 2009-10-02 06:43:51 kb Exp $
#
# This file contains rules which are shared between multiple Makefiles.
# This file is quite complicated - all compilation options are set in this
# file and included in all other Makefiles.
# Here you set the fortran compiler and the compilation mode (debug,prof,prod)
# If you want to include info for a new compiler - please know what you are 
# doing - then add a new file to the directory ../compilers.

SHELL   = /bin/sh
CPP	= /lib/cpp

# Top of this version of GETM and GOTM. Defaults can be overwritten using 
# environment variables.
ifndef GETMDIR
GETMDIR  = $(HOME)/GETM/getm-git
endif

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif

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

# OpenMP computation (mostly with compiler directives anyway):
ifeq ($(GETM_OMP),true)
DEFINES += -DGETM_OMP
endif

# Compile for parallel execution
ifeq ($(GETM_PARALLEL),true)
parallel=true
set par=par
else
parallel=false
set par=ser
endif

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

# Further directory related settings.
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

# Turbulence directory
GOTMLIBDIR	= $(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LINKDIRS	+= -L$(GOTMLIBDIR)
EXTRA_LIBS	+= -lturbulence$(buildtype) -lutil$(buildtype) 
INCDIRS		+= -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)

# Where does the NetCDF include file and library reside.

ifeq ($(NETCDF_VERSION),NETCDF4)

DEFINES		+= -DNETCDF4
INCDIRS		+= $(shell nc-config --fflags)
NETCDFLIB	=  $(shell nc-config --flibs)

else  # NetCDF3 is default

DEFINES		+= -DNETCDF3
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

endif

EXTRA_LIBS	+= $(NETCDFLIB)

# NetCDF/HDF configuration done

# Where does the MPI library reside.
ifeq ($(parallel),true)
DEFINES += -DGETM_PARALLEL

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

ifdef GETM_OMP
FLAGS   += $(OMP_FLAGS)
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
ifdef GETM_OMP
LDFLAGS += $(OMP_FLAGS)
endif

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
