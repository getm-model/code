#
# Makefile to build utilities written in Fortran90 - libfutils.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libfutils${buildtype}.a

MODSRC	= exceptions.F90 parallel.F90 parameters.F90 time.F90       \
          grid_interpol.F90

LIBSRC	= ver_interpol.F90 kbk_interpol.F90 tridiagonal.F90 pos.F90 \
          cnv_2d.F90 cnv_3d.F90 eta_mask.F90 col_interpol.F90 	    \
          to_2d_vel.F90 to_2d_u.F90 to_2d_v.F90                     \
	  to_3d_vel.F90 to_3d_uu.F90 to_3d_vv.F90                   \
	  to_velx.F90 to_vely.F90 tow.F90 c2x.F90 check_3d_fields.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD = ${LIB}(getm_timers.o)
ifeq ($(GETM_PARALLEL),true)
MOD += \
${LIB}(halo_mpi.o)	\
${LIB}(read_par_setup.o)
endif
MOD += \
${LIB}(exceptions.o)	\
${LIB}(parallel.o)	\
${LIB}(parameters.o)	\
${LIB}(halo_zones.o)	\
${LIB}(time.o)		\
${LIB}(getm_timers.o)	\
${LIB}(grid_interpol.o)

OBJ = \
${LIB}(ver_interpol.o)	\
${LIB}(kbk_interpol.o)	\
${LIB}(tridiagonal.o)	\
${LIB}(pos.o)		\
${LIB}(cnv_2d.o)	\
${LIB}(cnv_3d.o)	\
${LIB}(eta_mask.o)	\
${LIB}(col_interpol.o)	\
${LIB}(to_2d_vel.o)	\
${LIB}(to_2d_u.o)	\
${LIB}(to_2d_v.o)	\
${LIB}(to_3d_vel.o)	\
${LIB}(to_3d_uu.o)	\
${LIB}(to_3d_vv.o)	\
${LIB}(to_velx.o)	\
${LIB}(to_vely.o)	\
${LIB}(tow.o)		\
${LIB}(c2x.o)		\
${LIB}(check_3d_fields.o)	\
${LIB}(filter_1d.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/futils.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/{parameters,time,getm_timers,grid_interpol}.{m,mod}

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
