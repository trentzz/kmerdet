# ============================================================
# kmerdet paper — build with pdflatex + bibtex
# ============================================================

MAIN    = main
TEXFILE = $(MAIN).tex
PDFFILE = $(MAIN).pdf

# All chapter sources (order matches \input in main.tex)
CHAPTERS = chapters/introduction.tex \
           chapters/background.tex \
           chapters/methods.tex \
           chapters/implementation.tex \
           chapters/results.tex \
           chapters/discussion.tex \
           chapters/conclusion.tex

BIB = references.bib

.PHONY: all clean distclean view

# Default: build PDF
all: $(PDFFILE)

# Two-pass pdflatex + bibtex cycle required for citations and TOC
$(PDFFILE): $(TEXFILE) $(CHAPTERS) $(BIB)
	pdflatex -interaction=nonstopmode $(MAIN)
	bibtex $(MAIN)
	pdflatex -interaction=nonstopmode $(MAIN)
	pdflatex -interaction=nonstopmode $(MAIN)

# Remove intermediate build artifacts but keep the PDF
clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log \
	      $(MAIN).out $(MAIN).toc $(MAIN).lof $(MAIN).lot \
	      $(MAIN).loa $(MAIN).run.xml $(MAIN).bcf \
	      chapters/*.aux

# Remove everything including the PDF
distclean: clean
	rm -f $(PDFFILE)

# Open PDF with the system viewer
view: $(PDFFILE)
	xdg-open $(PDFFILE) 2>/dev/null || open $(PDFFILE) 2>/dev/null || evince $(PDFFILE)
