#
# Makefile to build the 2D specific library - libm2d.a
#

include ../Rules.make

INCS	= static_2d.h dynamic_declarations_2d.h dynamic_allocations_2d.h
LIB	= $(LIBDIR)/lib2d${buildtype}.a

TEXSRC	= m2d.F90 variables_2d.F90 \
	  advection.F90 adv_split_u.F90 adv_split_v.F90 adv_arakawa_j7_2dh.F90 adv_upstream_2dh.F90 adv_fct_2dh.F90 adv_tvd_limiter.F90 \
	  bdy_2d.F90 \
	  bottom_friction.F90 uv_advect.F90 uv_diffusion.F90 uv_diff_2dh.F90 momentum.F90 sealevel.F90 depth_update.F90 \
	  residual.F90 cfl_check.F90

MOD   = \
${LIB}(variables_2d.o)		\
${LIB}(bdy_2d.o)		\
${LIB}(advection.o)		\
${LIB}(m2d.o)

OBJ   = \
${LIB}(cfl_check.o)	\
${LIB}(bottom_friction.o)	\
${LIB}(depth_update.o)		\
${LIB}(momentum.o)		\
${LIB}(sealevel.o)		\
${LIB}(uv_advect.o)		\
${LIB}(uv_diffusion.o)		\
${LIB}(residual.o)		\
${LIB}(update_2d_bdy.o)		\
${LIB}(adv_split_u.o)		\
${LIB}(adv_split_v.o)		\
${LIB}(adv_upstream_2dh.o)	\
${LIB}(adv_fct_2dh.o)		\
${LIB}(adv_arakawa_j7_2dh.o)	\
${LIB}(adv_tvd_limiter.o)	\
${LIB}(uv_diff_2dh.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

$(MOD): $(INCS)

$(OBJ): $(INCS)

doc: $(TEXSRC)
	$(PROTEX) $(TEXSRC) > $(DOCDIR)/2d.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/m2d.{m,mod}

realclean: clean
	$(RM) *.o  doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
