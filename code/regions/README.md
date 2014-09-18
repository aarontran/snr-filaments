README, code
============

Aaron Tran
(last updated: 2014 September 14)

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

`ds9proj2box.py`

Converts ds9 region file with projections into ds9 file with boxes instead
of projections.  Then, region files w/ boxes can be safely converted to CIAO
files (using `regparse.conv_reg_coords`)

`regparse.py`

Module with useful region parsing methods

Plotting
--------

`fplot.py` and `matplotlibrc`

For nice and reproducible plots.  Ideally publication ready, although
in reality you probably have to twiddle plots individually.

Most useful is that plots and font sizes are pre-set for one column figures.


