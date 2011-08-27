#
# Master Makefile for making the 'getm' executable.
#
include Rules.make

#$(GETMDIR)/src/halos		\

# Set the subdirectories with Makefiles
SUBDIRS = \
$(GETMDIR)/src/futils		\
$(GETMDIR)/src/domain		\
$(GETMDIR)/src/meteo		\
$(GETMDIR)/src/2d
ifneq ($(GETM_NO_3D),true)
SUBDIRS += \
$(GETMDIR)/src/3d
endif
SUBDIRS += \
$(GETMDIR)/src/output		\
$(GETMDIR)/src/ncdf		\
$(GETMDIR)/src/input		\
$(GETMDIR)/src/getm

CORE_LIBS	=	\
		-lgetm${buildtype} 	\
		-loutput${buildtype}	\
		-ldomain${buildtype} 	\
		-lmeteo${buildtype}
ifneq ($(GETM_NO_3D),true)
CORE_LIBS	+=	\
		-l3d${buildtype}
endif
CORE_LIBS	+=	\
		-l2d${buildtype}	\
		-ldomain${buildtype} 	\
		-linput${buildtype}	\
		-lncdfio${buildtype}	\
		-lfutils${buildtype}

#		-lhalo${buildtype} 	\

LIBS		= $(CORE_LIBS) $(EXTRA_LIBS)

all: ../VERSION getm$(buildtype)_$(FORTRAN_COMPILER) install

../VERSION: ../Makefile Makefile Rules.make
	$(MAKE) -C ..

getm$(buildtype)_$(FORTRAN_COMPILER): makedirs subdirs
	$(FC) -o $@ getm/main.o $(LDFLAGS) $(LIBS)
	$(RM) getm/main.o

install: getm$(buildtype)_$(FORTRAN_COMPILER)
	mkdir -p $(BINDIR)
	mv getm$(buildtype)_$(FORTRAN_COMPILER) $(BINDIR)

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
	$(RM) ../bin/getm_*

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
