#$Id: Makefile,v 1.1.1.1 2002-05-02 14:01:56 gotm Exp $
# 
# Makefile to build the turbulence library - libturb.a.
#

include ../Rules.make

LIB	= $(LIBDIR)/libturb${buildtype}.a

MODSRC 	= turb.F90

LIBSRC 	= dissipation.F90 pran_kol.F90 production.F90 tke_calc.F90 eddyviscosity.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD	= \
${LIB}(turb.o)

OBJ	= \
${LIB}(dissipation.o)	\
${LIB}(pran_kol.o)	\
${LIB}(production.o)	\
${LIB}(tke_calc.o)	\
${LIB}(eddyviscosity.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/turb.tex

clean:
	$(RM) $(LIB) $(MODDIR)/turb.{m,mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
