#$Id: Makefile,v 1.1 2002-05-02 14:01:11 gotm Exp $
# 
# Makefile to build the domain library - libdomain.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libdomain${buildtype}.a

MODSRC	= domain.F90

LIBSRC	= bdy_spec.F90 print_bdy.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(domain.o)

OBJ   = \
${LIB}(bdy_spec.o)	\
${LIB}(print_bdy.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/domain.tex

clean:
	$(RM) $(LIB) $(MODDIR)/domain.{m.mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
