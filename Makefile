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

TAG=v$(shell echo $VER | tr . _)
RELEASE=getm-$(VER)
release_dir=/public/ftp/pub/getm-releases

RHOST=bolding-burchard.com
RUSER=bbh
RDIR=.

RHOST=gate
RUSER=kbk
RDIR=bolding-burchard.com/

base =  BUGS INSTALL README TODO HISTORY

.PHONY: doc

all: VERSION

VERSION: Makefile
	$(MAKE) distclean
	@echo $(VER) > $@
	@date > timestamp
	@echo \#define RELEASE \"$(VER)\" > .ver
	@mv -f .ver include/version.h

Makefile:

doc:
	$(MAKE) -C doc/

clean:
	rm -f VERSION

distclean:
	$(MAKE) -C src distclean
#	$(MAKE) -C utils distclean
	$(RM) -r bin/

tag:
	cvs tag $(TAG)

devel: tag
	(cd $(release_dir)/$@/ ; cvs -d gate:/public/cvs export -r $(TAG) -d getm-$(VER)/ getm-src )
	cvs2cl -F trunk
	mv ChangeLog $(release_dir)/$@/getm-$(VER)/
	(cd $(release_dir)/$@ ; tar -cvzf getm-$(VER).tar.gz getm-$(VER)/ )
	sync
	scp $(release_dir)/$@/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):$(RDIR)/src/$@/ChangeLog
	scp $(release_dir)/$@/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):$(RDIR)/src/$@/
	ssh $(RUSER)@$(RHOST) \( cd $(RDIR)/src/$@ \; \
	     ln -sf getm-$(VER).tar.gz getm-devel.tar.gz \)

stable: export
	(cd $(release_dir)/$@/ ; cvs -d gate:/public/cvs export -r $(TAG) -d getm-$(VER)/ getm-src )
	cvs2cl -b -F $(BRANCH) --no-ancestors
	mv ChangeLog $(release_dir)/$@/getm-$(VER)/
	(cd $(release_dir)/$@ ; tar -cvzf getm-$(VER).tar.gz getm-$(VER)/ )
	sync
	scp $(release_dir)/$@/getm-$(VER)/ChangeLog \
	    $(RUSER)@$(RHOST):src/$@/ChangeLog
	scp $(release_dir)/$@/getm-$(VER).tar.gz \
	    $(RUSER)@$(RHOST):src/$@/
	ssh $(RUSER)@$(RHOST) \( cd src/$@/ \; \
	     ln -sf getm-$(VER).tar.gz getm-stable.tar.gz \)

dummy:
