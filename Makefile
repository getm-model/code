#
# Master Makefile for the getm project.
#

VERSION = 1
PATCHLEVEL = 1
SUBLEVEL = 0

VER	= $(VERSION).$(PATCHLEVEL).$(SUBLEVEL)

EXEC	= model$(libtype)

base =  BUGS INSTALL README TODO HISTORY

all: VERSION $(EXEC) doc

$(EXEC): include/version.h
	$(MAKE) -C src $(EXEC) install

include/version.h: ./Makefile
	@echo \#define RELEASE \"$(VERSION).$(PATCHLEVEL).$(SUBLEVEL)\" > .ver
	@mv -f .ver $@

VERSION: ./Makefile
	 @echo $(VERSION).$(PATCHLEVEL).$(SUBLEVEL) > $@

doc:
	$(MAKE) -C src/ doc
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i doc; done

clean:
	rm -f VERSION

distclean:
	$(MAKE) -C src distclean
	$(MAKE) -C utils distclean
	$(RM) -r bin

dist: distclean
	(cd ../ ; tar cf - v$(VER)/ | gzip -9 > v$(VER).tar.gz) 
	sync

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
