#$Id: Makefile,v 1.5 2005-04-25 07:55:49 kbk Exp $
# 
# Makefile to build utilities written in Fortran90 - libfutils.a 
#

include ../Rules.make

LIB	= $(LIBDIR)/libfutils${buildtype}.a

MODSRC	= exceptions.F90 parallel.F90 parameters.F90 time.F90 grid_interpol.F90

LIBSRC	= ver_interpol.F90 kbk_interpol.F90 tridiagonal.F90 pos.F90 \
          cnv_2d.F90 cnv_3d.F90 eta_mask.F90 col_interpol.F90 	    \
          to_2d_vel.F90 to_3d_vel.F90 tow.F90

SRC	= $(MODSRC) $(LIBSRC)

ifeq ($(parallel),false)
MOD =
else
MOD = ${LIB}(halo_mpi.o) ${LIB}(read_par_setup.o)
endif

MOD += \
${LIB}(exceptions.o)	\
${LIB}(parallel.o)	\
${LIB}(parameters.o)	\
${LIB}(halo_zones.o)	\
${LIB}(time.o)		\
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
${LIB}(to_3d_vel.o)	\
${LIB}(tow.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/futils.tex

clean:
	$(RM) $(LIB) $(MODDIR)/{parameters,time,grid_interpol}.{m,mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
