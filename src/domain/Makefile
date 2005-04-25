#$Id: Makefile,v 1.3 2005-04-25 09:32:34 kbk Exp $
# 
# Makefile to build the domain library - libdomain.a
#

include ../Rules.make

INCS	= static_domain.h dynamic_declarations_domain.h 	\
          dynamic_allocations_domain.h
LIB	= $(LIBDIR)/libdomain${buildtype}.a

MODSRC	= topo_interface.F90 domain.F90 

LIBSRC	= part_domain.F90 bdy_spec.F90 mirror_bdy_2d.F90 print_bdy.F90 

SRC	= $(MODSRC) $(LIBSRC)

MOD   = \
${LIB}(topo_interface.o) \
${LIB}(domain.o)        


OBJ   = \
${LIB}(part_domain.o)	\
${LIB}(bdy_spec.o)	\
${LIB}(print_bdy.o)	\
${LIB}(mirror_bdy_2d.o)	\
${LIB}(mirror_bdy_3d.o)

all: modules objects

modules: $(MOD)

objects: $(OBJ)

$(MOD):	$(INCS)

$(OBJ):	$(INCS)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/domain.tex

clean:
	$(RM) $(LIB) $(MODDIR)/domain.{m.mod} $(MODDIR)/topo_interface.{m.mod}

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
