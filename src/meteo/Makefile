#$Id: Makefile,v 1.1 2002-05-02 14:01:38 gotm Exp $
# 
# Makefile to build the meteo module - meteo.mod
#

include ../Rules.make

MODSRC	= meteo.F90

LIBSRC	= exchange_coefficients.F90 short_wave_radiation.F90 fluxes.F90

SRC	= $(MODSRC) $(LIBSRC)

LIB	= $(LIBDIR)/libmeteo${buildtype}.a

MOD   = \
${LIB}(meteo.o)

OBJ   = \
${LIB}(short_wave_radiation.o)		\
${LIB}(exchange_coefficients.o)		\
${LIB}(fluxes.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/meteo.tex

clean:
	$(RM) $(LIB)

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
