#$Id: Makefile,v 1.2 2002-05-29 13:37:49 gotm Exp $
#
# Makefile to build the 3D specific library - libm3d.a
#

include ../Rules.make

INCS	= static_3d.h dynamic_declarations_3d.h dynamic_allocations_3d.h
LIB	= $(LIBDIR)/lib3d${buildtype}.a

MODSRC	= m3d.F90 variables_3d.F90 advection_3d.F90 eqstate.F90 \
          temperature.F90 salinity.F90 spm.F90 bdy_3d.F90 rivers.F90

LIBSRC	= start_macro.F90 bdy_3d.F90 coordinates.F90 bottom_friction_3d.F90 internal_pressure.F90 uu_momentum_3d.F90 vv_momentum_3d.F90 ww_momentum_3d.F90 uv_advect_3d.F90 uv_diffusion_3d.F90 stresses_3d.F90 ss_nn.F90 gotm.F90 slow_bottom_friction.F90 slow_advection.F90 slow_diffusion.F90 slow_terms.F90 stop_macro.F90

SRC	= $(MODSRC) $(LIBSRC)

MOD	=	\
$(LIB)(variables_3d.o)		\
$(LIB)(advection_3d.o)		\
$(LIB)(temperature.o)		\
$(LIB)(salinity.o)		\
$(LIB)(eqstate.o)		\
$(LIB)(spm.o)			\
$(LIB)(bdy_3d.o)		\
$(LIB)(m3d.o)			\
$(LIB)(rivers.o)

OBJ	=	\
$(LIB)(start_macro.o)		\
$(LIB)(coordinates.o)		\
$(LIB)(bottom_friction_3d.o)	\
$(LIB)(uu_momentum_3d.o)	\
$(LIB)(vv_momentum_3d.o)	\
$(LIB)(ww_momentum_3d.o)	\
$(LIB)(uv_advect_3d.o)		\
$(LIB)(uv_diffusion_3d.o)	\
$(LIB)(internal_pressure.o)	\
$(LIB)(eddyviscosity.o)		\
$(LIB)(slow_bottom_friction.o)	\
$(LIB)(slow_advection.o)	\
$(LIB)(slow_diffusion.o)	\
$(LIB)(slow_terms.o)		\
$(LIB)(stop_macro.o)

ifeq ($(turbulence),gotm)
OBJ += 		\
$(LIB)(stresses_3d.o)		\
$(LIB)(ss_nn.o)			\
$(LIB)(gotm.o)
endif

all: modules objects

modules: $(MOD)

objects: $(OBJ)

$(MOD): $(INCS)

$(OBJ): $(INCS)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/3d.tex

clean:
	$(RM) $(LIB) $(MODDIR)/m3d.{m.mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
