#$Id: Makefile,v 1.7 2007-09-21 13:03:42 kbk Exp $
#
# Makefile to build the ncdf module - ncdfio.o and the libncdfio.a library
#

include ../Rules.make

LIB	= $(LIBDIR)/libncdfio${buildtype}.a

MODSRC	= ncdf_common.F90 grid_ncdf.F90 ncdf_2d.F90 ncdf_2d_bdy.F90	\
          ncdf_3d.F90 ncdf_3d_bdy.F90 ncdf_meteo.F90 ncdf_rivers.F90	\
          ncdf_restart.F90 ncdf_mean.F90 ncdf_topo.F90

LIBSRC	= set_attributes.F90 init_grid_ncdf.F90 save_grid_ncdf.F90	\
          init_2d_ncdf.F90 save_2d_ncdf.F90				\
          init_3d_ncdf.F90 save_3d_ncdf.F90				\
          get_field_ncdf.F90 read_field_ncdf.F90 ncdf_close.F90		\
	  create_restart_ncdf.F90 write_restart_ncdf.F90		\
	  open_restart_ncdf.F90 read_restart_ncdf.F90			\
	  init_mean_ncdf.F90 save_mean_ncdf.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD	= 	                \
${LIB}(ncdf_common.o)		\
${LIB}(grid_ncdf.o)		\
${LIB}(ncdf_2d.o)		\
${LIB}(ncdf_2d_bdy.o)
ifneq ($(GETM_NO_3D),true)
MOD	+=                      \
${LIB}(ncdf_3d.o)		\
${LIB}(ncdf_3d_bdy.o)		\
${LIB}(ncdf_rivers.o)
endif
MOD	+= 	\
${LIB}(ncdf_meteo.o)		\
${LIB}(ncdf_topo.o)		\
${LIB}(ncdf_restart.o)		\
${LIB}(ncdf_mean.o)		

OBJ	= 	\
${LIB}(set_attributes.o)	\
${LIB}(init_grid_ncdf.o)	\
${LIB}(save_grid_ncdf.o)	\
${LIB}(init_2d_ncdf.o)		\
${LIB}(save_2d_ncdf.o)
ifneq ($(GETM_NO_3D),true)
OBJ	+= 	\
${LIB}(init_3d_ncdf.o)		\
${LIB}(save_3d_ncdf.o)		\
${LIB}(create_restart_ncdf.o)	\
${LIB}(write_restart_ncdf.o)	\
${LIB}(open_restart_ncdf.o)	\
${LIB}(read_restart_ncdf.o)	\
${LIB}(init_mean_ncdf.o)	\
${LIB}(save_mean_ncdf.o)	\
${LIB}(get_field_ncdf.o)	\
${LIB}(read_field_ncdf.o)
endif
OBJ	+= 	\
${LIB}(ncdf_close.o)

all: modules objects 

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/ncdf.tex
	touch doc

clean:
	$(RM) $(LIB) 

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
