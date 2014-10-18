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
${LIB}(cleanup.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/getm.tex
	touch doc

clean:
	$(RM) ${LIB} $(MODDIR)/{integration,initialise}.{m.mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean
	$(RM) work.pc *.d

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
