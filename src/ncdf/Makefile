#$Id: Makefile,v 1.3 2004-03-29 15:46:08 kbk Exp $
#
# Makefile to build the ncdf module - ncdfio.o and the libncdfio.a library
#

include ../Rules.make

LIB	= $(LIBDIR)/libncdfio${buildtype}.a

MODSRC	= ncdf_in.F90 ncdf_common.F90 ncdf_2d.F90 ncdf_2d_bdy.F90 \
          ncdf_3d.F90 ncdf_3d_bdy.F90 ncdf_meteo.F90 ncdf_rivers.F90 \
          ncdf_mean.F90 

LIBSRC	= ncdf_topo.F90 set_attributes.F90 init_2d_ncdf.F90 save_2d_ncdf.F90 \
      init_3d_ncdf.F90 save_3d_ncdf.F90 get_field_ncdf.F90 read_field_ncdf.F90 \
       ncdf_close.F90 init_mean_ncdf.F90 save_mean_ncdf.F90 towas.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD	= 	\
${LIB}(ncdf_common.o)		\
${LIB}(ncdf_2d.o)		\
${LIB}(ncdf_2d_bdy.o)
ifneq ($(GETM_NO_3D),true)
MOD	+= 	\
${LIB}(ncdf_3d.o)		\
${LIB}(ncdf_3d_bdy.o)		\
${LIB}(ncdf_rivers.o)
endif
MOD	+= 	\
${LIB}(ncdf_meteo.o)		\
${LIB}(ncdf_mean.o)		\
${LIB}(ncdf_in.o)

OBJ	= 	\
${LIB}(ncdf_topo.o)		\
${LIB}(set_attributes.o)	\
${LIB}(init_2d_ncdf.o)		\
${LIB}(save_2d_ncdf.o)
ifneq ($(GETM_NO_3D),true)
OBJ	+= 	\
${LIB}(init_3d_ncdf.o)		\
${LIB}(save_3d_ncdf.o)		\
${LIB}(init_mean_ncdf.o)	\
${LIB}(save_mean_ncdf.o)	\
${LIB}(get_field_ncdf.o)	\
${LIB}(towas.o)			\
${LIB}(read_field_ncdf.o)
endif
OBJ	+= 	\
${LIB}(ncdf_close.o)

all: modules objects 

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/ncdf.tex

clean:
	$(RM) $(LIB) $(MODDIR)/{ncdf_in,ncdf_out}.{m,mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
