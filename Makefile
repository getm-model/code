#
# Master Makefile for the getm project.
#

# 2003/09/07
VER=1.1.0
# 2003/09/16
VER=1.1.1
# 2003/10/18
VER=1.1.2
# 2004/01/02
VER=1.1.3
# 2004/01/08
VER=1.1.4
# 2004/04/06
VER=1.1.5
# 2004/06/15
VER=1.1.6
# 2005/01/14
VER=1.1.7
# 2005/04/20
VER=1.2.0_branch
# 2005/04/20
VER=1.3.0
# 2005/04/27
VER=1.3.1
# 2005/05/25
VER=1.3.2
# 2006/03/10
VER=1.3.3
# 2006/03/24 - new stable branch
VER=1.4
# 2006/03/24
VER=1.5.0
# 2006/03/24 - branch for developing adaptive coordinates
VER=adaptive_cord
# 2006/08/25
VER=1.5.1
# 2007/05/14 - new stable release
VER=1.6.0
# 2007/05/14 - new devel release
VER=1.7.0
# 2010/04/01 - new stable release
VER=1.8.0
# 2010/04/01 - new devel release
VER=1.9.0
# 2011/04/01 - new stable release
VER=2.0.0
# 2011/04/01 - new devel release
VER=2.1.0
# 2012/04/01 - new stable release
VER=2.2.0
# 2012/04/01 - new devel release
VER=2.3.0
# 2012/06/26 - unified advection/diffusion
VER=2.3.1
# 2013/04/01 - new stable release
VER=2.4.0
# 2013/04/01 - new devel release
VER=2.5.0

include compilers/compiler.$(FORTRAN_COMPILER)

.PHONY: doc

all: VERSION

VERSION: Makefile src/Makefile src/Rules.make
	$(MAKE) distclean
	@echo $(VER) > $@
	@date > timestamp
	@echo \#define RELEASE \"$(VER)\" > .ver
	@mv -f .ver include/version.h

FORTRAN:
	@echo "#define FORTRAN_VERSION \"`$(FC) $(VERSION_FLAG) 2> /dev/null`\"" > ./include/fortran_version.h

GIT:
	@echo "#define GIT_REVISION \"`git log | head -1`\"" > ./include/git_revision.h

Makefile:

doc:
	$(MAKE) -C doc/

devel stable branch: VERSION
	@echo
	@echo "making a new "$@" release: v"$(VER)
	@echo
	./release.sh $@ $(VER)

clean:
	rm -f VERSION

distclean:
	$(MAKE) -C doc $@
	$(MAKE) -C src $@
	$(RM) timestamp VERSION include/version.h
	$(RM) -r bin/ lib/ modules/

#-----------------------------------------------------------------------
# Copyright (C) 2006 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
