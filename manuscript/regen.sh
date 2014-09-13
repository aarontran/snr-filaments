#!/bin/sh

pdflatex paper-tycho.tex
bibtex paper-tycho
pdflatex paper-tycho.tex
pdflatex paper-tycho.tex
open paper-tycho.pdf
