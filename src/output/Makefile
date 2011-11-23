# 
# Makefile to build the output module
#

include ../Rules.make

LIB	= $(LIBDIR)/liboutput${buildtype}.a

MODSRC	= output.F90 parser.F90 variable_info.F90 \
          ncdf_out.F90 ascii_out.F90 \
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

#$(LIB)(variable_info.o)		\
#$(LIB)(parser.o)		\

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

tests: test_parser test_varinfo

test_parser: modules objects test_parser.o
	$(FC) -o $@ $@.o $(LDFLAGS) $(LIB)
	$(RM) $@.o
	./$@

test_varinfo: modules objects test_varinfo.o
	$(FC) -o $@ $@.o $(LDFLAGS) $(LIB)
	$(RM) $@.o
	./$@

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
