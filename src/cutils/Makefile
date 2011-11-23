#
# Makefile to build utilities written in C - libmcutils.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libcutils${buildtype}.a

OBJ   = \
${LIB}(julday.o)		\
${LIB}(caldat.o)

all: $(OBJ)

doc: $(SRC)
	$(PROTEX) $(SRC) > $(DOCDIR)/cutils.tex
	touch doc

clean:
	$(RM) $(LIB)

realclean: clean
	$(RM) *.o doc

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
