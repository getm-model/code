#$Id: Makefile,v 1.1 2002-05-02 14:00:35 gotm Exp $
#
# Master Makefile for making the 'getm' executable.
#
include Rules.make

# Set the subdirectories with Makefiles
SUBDIRS = \
$(GETMDIR)/src/futils		\
$(GETMDIR)/src/halos		\
$(GETMDIR)/src/domain		\
$(GETMDIR)/src/meteo		\
$(GETMDIR)/src/2d		\
$(GETMDIR)/src/3d		\
$(GETMDIR)/src/output		\
$(GETMDIR)/src/ncdf		\
$(GETMDIR)/src/input		\
$(GETMDIR)/src/getm

ifneq ($(turbulence),gotm)
SUBDIRS += $(GETMDIR)/src/turb
endif


CORE_LIBS	=	\
		-lgetm${buildtype} 	\
		-lhalo${buildtype} 	\
		-loutput${buildtype}	\
		-ldomain${buildtype} 	\
		-lmeteo${buildtype}	\
		-l2d${buildtype}	\
		-l3d${buildtype}	\
		-ldomain${buildtype} 	\
		-linput${buildtype}	\
		-lncdfio${buildtype}	\
		-lfutils${buildtype}

LIBS		= $(CORE_LIBS) $(EXTRA_LIBS)

all: model$(buildtype)_$(FORTRAN_COMPILER) install

model$(buildtype)_$(FORTRAN_COMPILER): makedirs subdirs
	$(FC) -o $@ $(LDFLAGS) $(LIBS)

install: model$(buildtype)_$(FORTRAN_COMPILER)
	mv model$(buildtype)_$(FORTRAN_COMPILER) $(BINDIR)

.PHONY: subdirs $(SUBDIRS)

makedirs:
	mkdir -p $(MODDIR) $(LIBDIR)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

modules:
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

objects:
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

doc:
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done
	
bkup: distclean
	tar -czvf ../../backups/getm.$$USER.`date +%Y%m%d%H`_`cat ../VERSION`.tar.gz *

clean:

realclean: clean
	$(RM) ../modules/*.mod ../lib/lib*$(buildtype).a
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

distclean: realclean
	$(RM) ../lib/*/lib*.a
	$(RM) ../modules/*/*
	$(RM) ../bin/model_*

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
