#$Id: Makefile,v 1.2 2003-04-07 16:45:07 kbk Exp $
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
	$(RM) work.pc *.d 

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
