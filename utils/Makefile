#$Id: Makefile,v 1.1 2003-04-03 06:11:26 gotm Exp $
#---------------------------------------------------------------------+
#
#	Makefile for various getm utilities
#
#---------------------------------------------------------------------+

include $(GETMDIR)/src/Rules.make

LDFLAGS	= -L/opt/netcdf/lib -L$(GETMDIR)/lib/$(FORTRAN_COMPILER)
FLIBS	= -lfutils_prod $(NETCDFLIBNAME) -lm -lPEPCF90

FFLAGS	= -O -I$(NETCDFINC) -I$(GETMDIR)/modules/$(FORTRAN_COMPILER)
F90FLAGS= $(FFLAGS)

FOBJS	= tide.o gen_bdy_ncdf.o

PROGS	= make_tides

BINDIR	= ../bin

all: $(PROGS) install

make_tides: make_tides.o $(FOBJS)
	@ echo "Making $@ ..."
	$(FC) -o $@  $@.o $(LDFLAGS) $(FOBJS) $(FLIBS)

install: $(PROGS)
	mkdir -p $(BINDIR)
	mv $(PROGS) $(BINDIR)

doc: html man

html: $(PROG).pod 
	pod2html $(PROG).pod 
	$(MV) $(PROG).html $(htmldir)

man: $(PROG).pod 
	pod2man $(PROG).pod > $(PROG).1
	$(MV) $(PROG).1 $(mandir)

clean:
	$(RM) *.o

realclean: clean
	$(RM) -f  $(PROGS)

distclean: realclean
	$(RM) -f *.c~ *.h~

%.o: %.F
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@

