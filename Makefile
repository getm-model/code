#
# Master Makefile for the getm project.
#

# 2003/09/07
VER=1.1.0
# 2003/09/16
VER=1.1.1

TAG=v$(shell cat VERSION | tr . _)
RELEASE=getm-$(VERSION)

RHOST=bolding-burchard.com
RUSER=bbh
RDIR=.

RHOST=gate
RUSER=kbk
RDIR=/public/bolding-burchard.com/

EXEC	= model$(libtype)

base =  BUGS INSTALL README TODO HISTORY

all: VERSION $(EXEC) doc
	$(MAKE) -C src $(EXEC) install

include/version.h: ./Makefile
	@echo \#define RELEASE \"$(VER)\" > .ver
	@mv -f .ver $@

VERSION: ./Makefile
	 @echo $(VER) > $@

doc:
	$(MAKE) -C src/ doc
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i doc; done

clean:
	rm -f VERSION

distclean:
	$(MAKE) -C src distclean
	$(MAKE) -C utils distclean
	$(RM) -r bin/

export:
	(cd ~/getm-releases ; cvs export -r $(TAG) getm ; mv getm getm-$(VER)/)
	cvs2cl
	mv ChangeLog ~/getm-releases/getm-$(VER)/
	(cd ~/getm-releases ; tar -cvzf getm-$(VER).tar.gz getm-$(VER)/ )
	sync

devel: export
	scp ~/getm-releases/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):$(RDIR)/src/ChangeLog.devel
	scp ~/getm-releases/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):$(RDIR)/src/
	ssh $(RUSER)@$(RHOST) \( cd $(RDIR)/src \; \
	     ln -sf getm-$(VER).tar.gz getm_devel.tar.gz \)

stable:
	scp ~/getm-releases/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):src/ChangeLog
	scp ~/getm-releases/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):src/
	ssh $(RUSER)@$(RHOST) \( cd src \; \
	     ln -sf getm-$(VER).tar.gz getm_stable.tar.gz \)

diff:
	cvs diff > cvs.diff
	vi cvs.diff
	$(RM) cvs.diff

update:
	cvs update > cvs.update
	vi cvs.update
	$(RM) cvs.update
	cvs2cl
	vi ChangeLog

commit:
	cvs commit > cvs.commit
	vi cvs.commit

dummy:
