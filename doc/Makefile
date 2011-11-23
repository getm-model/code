#
# Makefile for the GETM source code documentation.
#

SHELL	   = /bin/sh

SRCDIR    = ../src

RUSER=bolding
RHOST=bolding-burchard.com
RDIR=/data/kamel/domains/bolding-burchard.com/doc/v1.8/


PROTEXSRC  = getm.tex		\
             domain.tex		\
             2d.tex		\
             3d.tex		\
             meteo.tex		\
             input.tex		\
             output.tex		\
             ncdf.tex		\
             futils.tex

TEXSRC     = introduction.tex	\
             getmIntro.tex	\
             domainIntro.tex	\
             2dIntro.tex	\
             3dIntro.tex	\
             meteoIntro.tex	\
             inputIntro.tex	\
             outputIntro.tex	\
             ncdfIntro.tex	\
             futilsIntro.tex


TEXHTMLSRC = introduction.tex      \
             mainIntro.tex         \
             meanflowIntroHTML.tex \
             turbulenceIntro.tex   \
             airseaIntro.tex       \
             observationsIntro.tex \
             outputIntro.tex       \
             utilIntroHTML.tex     \
             extraIntro.tex        \
             cases.tex             \
             gotm.bib

.PHONY: html clean realclean distclean upload $(PROTEXSRC)

all: a4.dvi html

a4: a4.dvi
	dvips a4.dvi
	dvipdf a4.dvi

a4.dvi: a4.tex $(TEXSRC) $(PROTEXSRC)
	latex  a4.tex
	bibtex a4
	latex  a4.tex
	latex  a4.tex
	latex  a4.tex

html: a4.tex $(TEXSRC) $(PROTEXSRC)
	mkdir -p html
	latex2html -dir html a4.tex

# Generate the protex dependend latex files (from the source code)
$(PROTEXSRC):
	$(MAKE) -C $(SRCDIR)/$(shell echo $@ | cut -f 1 -d .) doc

upload: a4 html
	scp -r a4.pdf $(RUSER)@$(RHOST):$(RDIR)/pdf/
	scp -r html/* $(RUSER)@$(RHOST):$(RDIR)/html/

clean:
	$(RM) *.log *.toc *.aux *.bbl *.blg *.bac *~

realclean: clean
	$(RM) $(PROTEXSRC)
	$(RM) *.dvi *.ps *.pdf
	$(RM) -r html/
ifdef PROTEXSRC
	set -e; for i in $(PROTEXSRC); do $(MAKE) -C $(SRCDIR)/$(shell echo $(i) | cut -f 1 -d .) $@; done
endif

distclean: realclean
	$(RM) *.gz

