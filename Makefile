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

TAG=v$(shell cat VERSION | tr . _)
RELEASE=getm-$(VERSION)

RHOST=bolding-burchard.com
RUSER=bbh
RDIR=.

RHOST=gate
RUSER=kbk
RDIR=bolding-burchard.com/

EXEC	= model$(libtype)

base =  BUGS INSTALL README TODO HISTORY

all: VERSION $(EXEC) doc
	$(MAKE) -C src $(EXEC) install

include/version.h: ./Makefile
	@echo \#define RELEASE \"$(VER)\" > .ver
	@mv -f .ver $@

VERSION: include/version.h
	 @echo $(VER) > $@

doc:
	$(MAKE) -C src/ doc
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i doc; done

clean:
	rm -f VERSION

distclean:
	$(MAKE) -C src distclean
#	$(MAKE) -C utils distclean
	$(RM) -r bin/

tag:
	cvs tag $(TAG)

export: tag
	(cd ~/getm-releases ; cvs export -r $(TAG) getm-src ; mv getm-src getm-$(VER)/)

devel: export
	cvs2cl -F trunk
	mv ChangeLog ~/getm-releases/getm-$(VER)/
	(cd ~/getm-releases ; tar -cvzf getm-$(VER).tar.gz getm-$(VER)/ )
	sync
	scp ~/getm-releases/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):$(RDIR)/src/ChangeLog.devel
	scp ~/getm-releases/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):$(RDIR)/src/
	ssh $(RUSER)@$(RHOST) \( cd $(RDIR)/src \; \
	     ln -sf getm-$(VER).tar.gz getm_devel.tar.gz \)

stable: export
	cvs2cl -b -F $(BRANCH) --no-ancestors
	mv ChangeLog ~/getm-releases/getm-$(VER)/
	(cd ~/getm-releases ; tar -cvzf getm-$(VER).tar.gz getm-$(VER)/ )
	sync
	scp ~/getm-releases/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):src/ChangeLog
	scp ~/getm-releases/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):src/
	ssh $(RUSER)@$(RHOST) \( cd src \; \
	     ln -sf getm-$(VER).tar.gz getm_stable.tar.gz \)

dummy:
