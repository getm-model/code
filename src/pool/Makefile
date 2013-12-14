#$Id: Makefile,v 1.6 2006-02-10 22:29:12 hb Exp $
#
# Makefile to build the POOL specific library - libpool.a
#

include ../Rules.make

INCS	=
LIB	= $(LIBDIR)/libpool${buildtype}.a

MODSRC	= pool.F90

LIBSRC	= deformation_rates.F90 flux_center2interface.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(pool.o)

OBJ   = \
${LIB}(deformation_rates.o)	\
${LIB}(flux_center2interface.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/pool.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/les.{m,mod}

realclean: clean
	$(RM) *.o  doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2013 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------