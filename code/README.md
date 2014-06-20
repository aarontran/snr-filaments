README for snr-research scripts
===============================

Aaron Tran
(last updated: 2014 June 18)

Spectrum manipulation
---------------------

`spec_linkbg.py`

Connects spectra for a list of regions, with closest background spectra
Reads CIAO region files used to generate 1) regions of interest, 2) backgrounds
in order to determine closest backgrounds
(future extension: allow to parse more types of CIAOREG regions)

`spec_clearbg.py`

Clears the BACKFILE entries for a set of spectra (specify the filename stem, it
will update the headers for ungrouped and grouped spectra).
Use this if you mess up the order of arguments to `spec_linkbg.py` and end up
writing headers to the background files.

`spec_fitplot.py`

Generate best fits (parameters, plots) for set of spectra, using PyXspec
Logs with fit parameters are saved to same folder as spectra
Must run using 32 bit Python (`arch -i386 python`)!!!!!  E.g.,
    arch -i386 python spec_fitplot.py ... stuff

`spec_plot2pdf.sh`

Short bash script to convert ps plots to pdf plots, deleting the ps plots in
the process

Region manipulation
-------------------

`ciaoreg2spec.py`

Wrapper function to run CIAO's specextract on CIAOREG files.
NOT YET IMPLEMENTED, NOT REALLY IMPORTANT.  Easy enough to
execute specextract from command line.

`ds9projplotter.py`

Takes a ds9 region file containing projections, generates output plots of the
projections.  Uses RGB data from Tycho, specifically.

`ds9proj2box.py`

Converts ds9 region file with projections into CIAO file, turning projections
into CIAO rotboxes.
(ideally, can be extended to convert more DS9 regions to CIAO equivalents)
(currently, it will only convert projections, and anything else that DS9
already knows how to convert)


Plotting (from my homebrewed software utilities)
------------------------------------------------

`fplot.py` and `matplotlibrc`

For nice and reproducible plots.  Ideally publication ready, although
in reality you probably have to twiddel plots individually.


