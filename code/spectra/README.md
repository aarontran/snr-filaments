Scripts/code for generating and manipulating Chandra spectra
============================================================

(last updated: 2014 September 14)

Nina/Brian's shell scripts to make/merge spectra from multiple ObsIDs
---------------------------------------------------------------------

### `newregion.sh`
Open ds9 window (ds9&) + call heainit before running script as:

    ./newregion.sh [NAME]

[NAME] should not include the `.reg` (file must end with `.reg`)

The original region file should be in DS9 fk5 coordinates;
this script generates output CIAO physical coordinates for EACH ObsID.
Due to slightly different quirks of each observation, these coordinates are in
general NOT THE SAME -- hence this script.

### `newspectrum.sh`
Call CIAO before running script:

    ./newspectrum.sh [NAME]

### `mergespectra.sh`
Call heainit and CIAO before running script (yes, need both)

After generating combined spectra with .wrmf, .warf files, uses HEASOFT
utility `grppha` to bin spectra if needed.

### `src_orig/`
The original files as given to me by Nina


My own utilities to manipulate spectra (some may be outdated)
-------------------------------------------------------------

`spec_linkbg.py`  (must run CIAO first)

Connects spectra for a list of regions, with closest background spectra.
Regions must be boxes, cannot be projections -- should use same region files
specified for newregion/newspectrum/mergespectrum shell scripts.
(future extension: allow to parse more types of CIAOREG regions)

`spec_clearbg.py`

Clears the BACKFILE entries for a set of spectra (specify the filename stem, it
will update the headers for ungrouped and grouped spectra).
Use this if you mess up the order of arguments to `spec_linkbg.py` and end up
writing headers to the background files.

`spec_fit.py`  (must run heainit first)

Generate best fits (parameters, plots) for set of spectra, using PyXspec
Logs with fit parameters are saved to same folder as spectra
Must run using 32 bit Python (`arch -i386 python`)!!!!!  E.g.,
    arch -i386 python spec_fitplot.py ... stuff

`spec_plot2pdf.sh`

Short bash script to convert ps plots to pdf plots, deleting the ps plots in
the process


