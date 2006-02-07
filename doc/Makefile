#$Id: Makefile,v 1.1 2006-02-07 07:16:21 kbk Exp $
#
# Makefile for the GETM source code documentation. 
#

SHELL	   = /bin/sh

SRCDIR    = ../src

RUSER=bolding
RHOST=bolding-burchard.com
RDIR=/data/kamel/domains/bolding-burchard.com/doc/v1.6/pdf/


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

.PHONY :     clean realclean distclean upload $(PROTEXSRC)

all: a4.dvi a4.ps a4.pdf

a4.dvi: a4.tex $(TEXSRC) $(PROTEXSRC) figures
	latex  a4.tex
	bibtex a4
	latex  a4.tex
	latex  a4.tex
	latex  a4.tex

a4.ps: a4.dvi
	dvips a4.dvi

a4.pdf: a4.dvi
	dvipdf a4.dvi

# Generate the protex dependend latex files (from the source code)
$(PROTEXSRC):
	$(MAKE) -C $(SRCDIR)/$(shell echo $@ | cut -f 1 -d .) doc

upload: a4.pdf
	scp -r a4.pdf $(RUSER)@$(RHOST):$(RDIR)

figures:
	wget -nd -nH -q http://www.bolding-burchard.com/src/figures/getm-figs-1.6.x.gz
	tar -xvzf getm-figs-1.6.x.gz
	$(RM) getm-figs-1.6.x.gz

clean:
	$(RM) *.log *.toc *.aux *.bbl *.blg *.bac *~

realclean: clean
	$(RM) *.dvi *.ps *.pdf
	$(RM) $(PROTEXSRC)

distclean: realclean
	$(RM) *.gz
	$(RM) -r figures

