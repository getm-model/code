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

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/meteo.tex
	touch doc

clean:
	$(RM) $(LIB)

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
