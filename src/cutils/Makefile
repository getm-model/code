#$Id: Makefile,v 1.1.1.1 2002-05-02 14:01:04 gotm Exp $
# 
# Makefile to build utilities written in C - libmcutils.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libcutils${buildtype}.a

OBJ   = \
${LIB}(julday.o)		\
${LIB}(caldat.o)

all: $(OBJ)

doc:
	$(PROTEX) $(SRC) > $(DOCDIR)/cutils.tex

clean:
	$(RM) $(LIB)

realclean: clean
	$(RM) *.o

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
