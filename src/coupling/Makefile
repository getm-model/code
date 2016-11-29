#
# Makefile to build the getm_esmf library - libgetm_esmf.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libgetm_esmf${buildtype}.a

SRC	= getm_esmf.F90

MOD   = \
${LIB}(getm_esmf.o)

OBJ   =

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/getm_esmf.tex
	touch doc

clean:
	$(RM) ${LIB} $(MODDIR)/{getm_esmf}.{m.mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
