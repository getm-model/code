#$Id: Makefile,v 1.5 2009-04-28 06:49:11 kb Exp $
# 
# Makefile to build the output module
#

include ../Rules.make

LIB	= $(LIBDIR)/liboutput${buildtype}.a

MODSRC	= output.F90 ncdf_out.F90 ascii_out.F90 \
          diagnostic_variables.F90 nesting.F90 

LIBSRC	= diagnose.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(nesting.o)
ifneq ($(GETM_NO_3D),true)
MOD += \
${LIB}(diagnostic_variables.o)
endif
MOD += \
${LIB}(ncdf_out.o)		\
${LIB}(ascii_out.o)		\
${LIB}(output.o)		

OBJ     =
ifneq ($(GETM_NO_3D),true)
OBJ += \
$(LIB)(calc_mean_fields.o)
endif
OBJ += \
$(LIB)(diagnose.o) 	

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/output.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/{output,ncdf_out,ascii_out}.{m,mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
