#
# Makefile to build the ncdf module - ncdfio.o
#

include ../Rules.make

LIB	= $(LIBDIR)/libinput${buildtype}.a

MODSRC	= input.F90

LIBSRC	= init_2d_bdy.F90 get_2d_bdy.F90 init_3d_bdy.F90 get_3d_bdy.F90 \
          init_meteo_input.F90 get_meteo_data.F90 \
          init_river_input.F90 get_river_data.F90 \
          get_2d_field.F90 get_3d_field.F90 \
          read_profile.F90 read_topo_file.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(input.o)

OBJ   = \
${LIB}(init_2d_bdy.o)		\
${LIB}(get_2d_bdy.o)
ifneq ($(GETM_NO_3D),true)
OBJ   += \
${LIB}(init_3d_bdy.o)		\
${LIB}(get_3d_bdy.o)		\
${LIB}(init_river_input.o)	\
${LIB}(get_river_data.o)
endif
OBJ   += \
${LIB}(init_meteo_input.o)	\
${LIB}(get_meteo_data.o)	\
${LIB}(get_2d_field.o)		\
${LIB}(get_3d_field.o)		\
${LIB}(read_profile.o)          \
${LIB}(read_topo_file.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/input.tex
	touch doc

clean:
	$(RM) $(LIB)  $(MODDIR)/input.{m,mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
