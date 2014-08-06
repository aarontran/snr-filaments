README, code
============

Aaron Tran
(last updated: 2014 August 6)

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

`spec_fit.py`

Generate best fits (parameters, plots) for set of spectra, using PyXspec
Logs with fit parameters are saved to same folder as spectra
Must run using 32 bit Python (`arch -i386 python`)!!!!!  E.g.,
    arch -i386 python spec_fitplot.py ... stuff

`spec_plot2pdf.sh`

Short bash script to convert ps plots to pdf plots, deleting the ps plots in
the process

`ciaoreg2spec.py`

Wrapper function to run CIAO's specextract on CIAOREG files.
NOT YET IMPLEMENTED, NOT REALLY IMPORTANT.  Easy enough to
execute specextract from command line.


Region manipulation
-------------------

`ds9projplotter.py`

Takes a ds9 region file containing projections, generates one or both of
(1) plots of projection profiles, (2) plaintext projection data.
Operates on an arbitrary number of fits files (bands).  Labels identify output
files to corresponding bands (but must be manually input at this time).

`ds9projsplitter.py`

Goes through pickled region dictionary from profile fits and splits regions
into two pieces, depending on FWHM locations.  This gives two region files, one
with upstream boxes and one with downstream boxes, which can be used to
generate spectra.

`ds9proj2ciao.py`

Converts ds9 region file with projections into CIAO file with rotboxes instead
of projections.  Other regions are converted to CIAO regions automatically
(i.e., some will convert as intended by DS9, others will drop out).
(ideally, extend script to convert more DS9 regions to CIAO equivalents)

`regparse.py`

Module with useful region parsing methods

Plotting
--------

`fplot.py` and `matplotlibrc`

For nice and reproducible plots.  Ideally publication ready, although
in reality you probably have to twiddle plots individually.


