#$Id: Makefile,v 1.1 2002-05-02 14:01:25 gotm Exp $
#
# Makefile to build the main getm library - libmgetm.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libgetm${buildtype}.a

SRC	= main.F90 initialise.F90 integration.F90 cleanup.F90

MOD   = \
${LIB}(integration.o)		\
${LIB}(initialise.o)

OBJ   = \
${LIB}(cleanup.o)		\
${LIB}(main.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/getm.tex

clean:
	$(RM) ${LIB} $(MODDIR)/{integration,initialise}.{m.mod}

realclean: clean
	$(RM) *.o 

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
