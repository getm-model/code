#$Id: Rules.make,v 1.13.8.2 2007-05-14 12:37:53 kbk Exp $
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

# The compilation mode is obtained from $COMPILATION_MODE
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

DEFINES=-D$(FORTRAN_COMPILER)

ifeq ($(GETM_NO_3D),true)
DEFINES += -DNO_3D
export GETM_NO_BAROCLINIC=true
endif

ifeq ($(GETM_NO_BAROCLINIC),true)
DEFINES += -DNO_BAROCLINIC
endif

ifeq ($(GETM_SPM),true)
DEFINES += -DSPM
endif

ifeq ($(GETM_BIO),true)
DEFINES += -DGETM_BIO
endif

# The compilation mode is obtained from $COMPILATION_MODE
# default production - else debug or profiling
ifndef GETM_PARALLEL
parallel=false
set par=ser
else
parallel=true
set par=par
endif

turbulence=
turbulence=gotm

CPP	= /lib/cpp

# Here you can put defines for the [c|f]pp - some will also be set depending
# on compilation mode - if STATIC is defined be careful.
#DEFINES = -DHAIDVOGEL_TEST
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
#DEFINES += -ECMWF_FRV
#DEFINES += -HIRLAM_FRV
#DEFINES += -SETTING_LON_LAN
ifdef STATIC
else
endif

# Directory related settings.

# Top of this version of getm.
ifndef GETMDIR
GETMDIR  = $(HOME)/GETM/getm-stable
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
ifdef NETCDFLIBNAME
NETCDFLIB	= $(NETCDFLIBNAME)
else
NETCDFLIB	= -lnetcdf
ifdef NETCDFLIBDIR
LINKDIRS	+= -L$(NETCDFLIBDIR)
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
ifdef MPILIBFILE
MPILIB		= $(MPILIBFILE)
else
ifdef MPILIBNAME
MPILIB		= $(MPILIBNAME)
else
MPILIB		= -lmpich
endif
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

#KBKinclude $(GETMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

#DEFINES += -DREAL_4B=$(REAL_4B)

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

include $(GETMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

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

