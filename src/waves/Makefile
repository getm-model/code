#$Id: Makefile,v 1.6 2006-02-10 22:29:12 hb Exp $
#
# Makefile to build the WAVES specific library - libwaves.a
#

include ../Rules.make

INCS	=
LIB	= $(LIBDIR)/libwaves${buildtype}.a

MODSRC	= waves.F90 variables_waves.F90

LIBSRC	= 

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(variables_waves.o)	\
${LIB}(waves.o)

OBJ   = \
${LIB}(stokes_drift.o)		\
${LIB}(rs_force.o)		\
${LIB}(radiation_stress.o)	\
${LIB}(vortex_force.o)		\
${LIB}(stokes_drift_3d.o)	\
${LIB}(radiation_stress_3d.o)	\
${LIB}(vortex_force_3d.o)	\
${LIB}(bottom_friction_waves.o)


all: modules objects

modules: $(MOD)

objects: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/waves.tex
	touch doc

clean:
	$(RM) $(LIB) $(MODDIR)/les.{m,mod}

realclean: clean
	$(RM) *.o  doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2013 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------