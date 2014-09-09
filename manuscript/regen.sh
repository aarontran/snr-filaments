#!/bin/sh

pdflatex tycho-paper.tex
bibtex tycho-paper
pdflatex tycho-paper.tex
pdflatex tycho-paper.tex
open tycho-paper.pdf
