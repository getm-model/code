#$Id: Makefile,v 1.1.1.1 2002-05-02 14:01:52 gotm Exp $
# 
# Makefile to build the output module
#

include ../Rules.make

LIB	= $(LIBDIR)/liboutput${buildtype}.a

MODSRC	= output.F90 ncdf_out.F90 ascii_out.F90 nesting.F90

LIBSRC	= diagnose.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(nesting.o)		\
${LIB}(ncdf_out.o)		\
${LIB}(ascii_out.o)		\
${LIB}(output.o)

OBJ     =       \
$(LIB)(diagnose.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/output.tex

clean:
	$(RM) $(LIB) $(MODDIR)/{output,ncdf_out,ascii_out}.{m,mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
