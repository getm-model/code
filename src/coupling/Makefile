#
# Makefile to build the coupling libraries - libgetm_esmf.a libgetm_oasis.a
#

include ../Rules.make

LIB_ESMF	= $(LIBDIR)/libgetm_esmf${buildtype}.a
LIB_OASIS	= $(LIBDIR)/libgetm_oasis${buildtype}.a
LIB		= $(LIB_ESMF) $(LIB_OASIS)

SRC	=
ifeq ($(GETM_ESMF),true)
SRC	+= getm_esmf.F90
endif
ifeq ($(GETM_OASIS),true)
SRC	+= getm_oasis.F90
endif

MOD   =
ifeq ($(GETM_ESMF),true)
MOD   += ${LIB_ESMF}(getm_esmf.o)
endif
ifeq ($(GETM_OASIS),true)
MOD   += ${LIB_OASIS}(getm_oasis.o)
endif

OBJ   =

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/coupling.tex
	touch doc

clean:
	$(RM) ${LIB} $(MODDIR)/{getm_esmf,getm_oasis}.{m.mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
