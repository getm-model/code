#$Id: Makefile,v 1.1.1.1 2002-05-02 14:01:29 gotm Exp $
# 
# Makefile to build the mpi interface library - libhalo.a
#

include ../Rules.make

#DEFINES += -DSTATIC

LIB	= $(LIBDIR)/libhalo${buildtype}.a

MODSRC	= parallel.F90 halo_zones.F90 halo_mpi.F90

LIBSRC	=

SRC	= $(MODSRC) $(LIBSRC)

MOD	= \
${LIB}(parallel.o)	\
${LIB}(halo_mpi.o)	\
${LIB}(halo_zones.o)

OBJ	= \
${LIB}(halo_test.o)

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
