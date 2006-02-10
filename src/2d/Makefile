#$Id: Makefile,v 1.6 2006-02-10 22:29:12 hb Exp $
#
# Makefile to build the 2D specific library - libm2d.a
#

include ../Rules.make

INCS	= static_2d.h dynamic_declarations_2d.h dynamic_allocations_2d.h
LIB	= $(LIBDIR)/lib2d${buildtype}.a

MODSRC	= m2d.F90 variables_2d.F90 

LIBSRC	= momentum.F90 uv_advect.F90  uv_diffusion.F90 bottom_friction.F90 sealevel.F90 depth_update.F90 uv_depths.F90 update_2d_bdy.F90  cfl_check.F90 have_bdy.F90 residual.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(variables_2d.o)		\
${LIB}(m2d.o)

OBJ   = \
${LIB}(cfl_check.o)	\
${LIB}(bottom_friction.o)	\
${LIB}(depth_update.o)		\
${LIB}(momentum.o)		\
${LIB}(sealevel.o)		\
${LIB}(uv_advect.o)		\
${LIB}(uv_diffusion.o)		\
${LIB}(uv_depths.o)		\
${LIB}(have_bdy.o)		\
${LIB}(residual.o)		\
${LIB}(update_2d_bdy.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

$(MOD): $(INCS)

$(OBJ): $(INCS)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/2d.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/m2d.{m,mod}

realclean: clean
	$(RM) *.o  doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
