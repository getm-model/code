#$Id: Makefile,v 1.2 2002-05-02 15:04:42 gotm Exp $
# 
# Makefile to build the mpi interface library - libhalo.a
#

include ../Rules.make

#DEFINES += -DSTATIC

LIB	= $(LIBDIR)/libhalo${buildtype}.a

MODSRC	= parallel.F90 halo_zones.F90 halo_mpi.F90

LIBSRC	=

SRC	= $(MODSRC) $(LIBSRC)

ifeq ($(parallel),false)
MOD	= \
${LIB}(parallel.o)
else
MOD	= \
${LIB}(halo_mpi.o)	\
${LIB}(halo_zones.o)
endif

ifeq ($(parallel),false)
OBJ	=
else
OBJ	= \
${LIB}(halo_test.o)
endif

all: modules objects

modules: $(MOD)

objects: $(OBJ)
#objects:

halo_test: modules objects
	$(FC) -o $@ $(LIB) -L$(MPILIBDIR) -lmpich

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/halo.tex

clean:
	$(RM) $(LIB) $(MODDIR)/halo.{m.mod}
	$(RM) ??-{domain,west,north,east,south} halo_test

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
