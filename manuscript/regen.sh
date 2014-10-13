#!/bin/sh

# Run bibtex if ANY command line arg is given
if [ "$#" -eq 1 ]
then
    pdflatex paper-tycho.tex
    bibtex paper-tycho
    pdflatex paper-tycho.tex
fi

pdflatex paper-tycho.tex
open paper-tycho.pdf
