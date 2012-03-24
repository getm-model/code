#$Id: Makefile,v 1.6 2006-02-10 22:29:12 hb Exp $
#
# Makefile to build the LES specific library - libles.a
#

include ../Rules.make

INCS	=
LIB	= $(LIBDIR)/libles${buildtype}.a

MODSRC	= les.F90 variables_les.F90

LIBSRC	= les_smagorinsky.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(variables_les.o)	\
${LIB}(les.o)

OBJ   = \
${LIB}(les_smagorinsky.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/les.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/les.{m,mod}

realclean: clean
	$(RM) *.o  doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
