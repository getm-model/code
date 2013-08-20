#
# Makefile to build the GETM ice module - getm_ice
#

include ../Rules.make

MODSRC	= getm_ice.F90

SRC	= $(MODSRC) $(LIBSRC)

LIB	= $(LIBDIR)/libgetm_ice${buildtype}.a

MOD   = \
${LIB}(getm_ice.o)

all: modules

modules: $(MOD)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/getm_ice.tex
	touch doc

clean:
	$(RM) $(LIB)

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2013 - Hans Burchard and Karsten Bolding (BB)          !
#-----------------------------------------------------------------------
