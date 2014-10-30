snr-research code pipeline
==========================
Aaron Tran
Summer 2014

N.B. this may not be fully up to date -- see `rsch-notes.md`, and talk to me
before using pipeline.  The hope is that it's straightforward enough to make
sense of, at least (see code documentation / command-line help).

The most fuzzy/hand-wavy part is region selection, see notes in data/ for that.
A lot of tweaking goes into 1. FWHM fitting, and 2. spectrum fitting.  Look at
the code to see what numbers are being set/frozen/thawed/etc.

Below I give an abbreviated version of the pipeline, for quick reference and to
see where files go / what code does what.

WARNING: the scripts aren't super robust -- they make assumptions about
filenames, numbers/ordering, etc..., which generally hold if you follow the
pipeline.  Manually moving files a lot is not the best idea (best to rerun the
pipeline).  Numbering/ordering follows region file order, starting from 1.

WARNING 2: always specify labels in order (0.7-1kev, 1-2kev, 2-7kev, etc...)
Many bits of code will assume ordered labels as they pass data around.

Directory structure
===================
Idealized, update as needed... (set-up 2014 July 8, somewhat following
the idea from a Software Carpentry lecture
[link](http://software-carpentry.org/v4/data/mgmt.html).

Ordered in order of generation, roughly (to make data dependencies clearer).
Paper tables are copy-pasted from IPython notebooks (computed from these data).
Paper plots in `region-ID/fwhms/plots` and `region-ID/fwhms/model-fits/plots`.

    data-*/ *.fits
            README.md

            region-ID/ regionID.reg
                       imgs/ *.jpg  # Images of regions on SNR
                             *.eps

                       profiles/ prf_*.dat
                                 prf-cts_*.dat
                                 prf-proc_*.dat
                                 prf-proc_*.npz  # These are stable
                                 plots/ *.png

                       fwhms/ fwhms_*.pkl
                              fwhms_*.json
                              fwhms_spec_cut.npz
                              fwhms_spec_cut.dat

                              plots/ prfs_*.pdf  # Fitted profiles/spectra
                                     spec_*.pdf

                              model-fits/ *.pkl  # All model fit results!
                                          *.json
                                          plots/ *prfs*.pdf
                                                 *prfs-radio*.pdf
                                                 *energywidth*.pdf

                        fwhms-cap/ ...  # Same structure
                        fwhms-nobkg/ ...  # BUT, spectra must always be
                                          # derived from fwhms/*

                       spectra/ *.pi  # Both down/up sets
                                *.rmf
                                *.arf
                                plots/ *.pdf   # Plots of spectra/fits
                                       *.ps
                                fits/ *.txt  # Fitting information
                                      *.qdp  # from XSPEC, etc.

            background-ID/ *.reg
                           *.ciaoreg
                           spectra/ *.pi
                                    *.rmf
                                    *.arf

    code/ fplot.py
          matplotlibrc

          regions/ ds9projplotter.py
                   ds9projsplitter.py
                   ds9proj2box.py
                   ds9projangle.py
                   (regparse.py)
                   (fplot.py, matplotlibrc)

          profiles/ prep_profile_fit.py
                    profile_process_*.ipynb
                    regdict.py
                    (crvfit.py, ffwhm.py, fit_funcs.py, fsmooth.py)

          spectra/ newregion.sh
                   newspectrum.sh
                   mergespectra.sh
                   spec_linkbg.py
                   spec_fit.py
                   (spec_clearbg.py)  # For when you mess up

          models/ *-fit-tables.ipynb
                  models_disp.py  # Still being twiddled
                  snr_catalog.py  # Contains important information
                  tables/...  # SNR specific tables -- but this is dependent on
                              # FWHM measurements!
                  (models.py, models_exec.py, parutils.py, full model codes)

          prettify/

    # Parentheses denote utility/helper-type modules

The key outputs are profiles, fitted FWHMs, spectra (based on FWHM fitting),
averaged FWHMs, and model fits.


Step-by-step guide to pipeline
==============================

I have tried to highlight locations where the user might have to intervene,
twiddle numbers, edit code, etc.

Pipeline 1: radial profiles and fits
====================================

## Prerequisites
Use CIAO's `reproject_obs`, `flux_obs`, and `dmcopy` to make energy band
images, in both intensity flux units and in uncorrected counts
(no exposure/vignetting correction).

## Region creation, cataloging
Open RGB image of SNR data.

    ./ds9-tycho-rgb.sh

Select regions of interest and backgrounds by hand.  Add text labels
(they will not show up in the CIAO regions).  Color code regions.
Pick regions with slight thermal uptick in back.  Save copies in fk5 and in
physical coordinates.

Output: `regions-n.reg`, `regions-n.physreg`

## Making nice images of regions

Open RGB image of SNR data.  Twiddle colormap/scaling parameters as desired.
Convert projection regions to boxes (for image only) with:

    code/regions/ds9proj2box.py

Remove automatically placed numbers/labels and add your own, using Times font.
Save the regions as `regions-n-disp.reg` or whatever.
It may be good to save the scaling/colormap parameters too (just use ds9
command line syntax and write in by hand).

## Compute azimuth angles for regions (Tycho)

For shock velocity interpolation (if needed)
Pick some circle region to compute the azimuth angles (`data-SNR/circ/`)

Code:   `code/regions/ds9projangle.py`
Input:  `regions-n.physreg`, `circ.physreg`
Output: `data-SNR/regions-n/regions-n-az.txt`

## Generate, process radial profiles from regions, mosaics and count files

Code:   `code/regions/ds9projplotter.py`
Input:  `regions-n.reg`, `*_mosaic.fits`, `*_counts.fits`
Output: `profiles/prf_[...].dat`, `profiles/prf-cts_[...].dat`

Code:   `code/profiles/prep_profile_fit.py`
Input:  `regions-n.physreg`, `prf*.dat`
Output: `prf-proc*.[dat,npz]`, `prf-proc_fit_cuts.[dat,npz]`

Set smoothing window, binsize (known a priori), etc.

## Fit radial profiles and obtain FWHMs

Code:   `code/profiles/profile_process_*.ipynb`
Input:  `regions-n.reg`, `regions-n.physreg`, `profiles/prf-proc_*.dat`
Output: `fwhms/fwhms_*.[pkl,json]`, `fwhms/fwhms_spec_cuts.[npz,dat]`

Review code in IPython notebook and make any necessary changes, adjustments,
etc.

## Generate plots of profiles and fits for paper

Code:   `code/profiles/plotter-prfs-specs.ipynb`
Input:  `fwhms/fwhms_*.[pkl,json]`, `fwhms/fwhms_spec_cuts.[npz,dat]`
Output: `fwhms/plots/prfs_{:02d}.pdf`

Plots of profiles and data in all energy bands, optimized for manuscript
Create the directory `plots` before running.



Pipeline 2: extract and fit spectra
===================================
Prerequisite: 1

## Prep/split regions based on profile fits

This code handles the projection-box conversion as well

Code:   `../code/ds9projsplitter.py`
Input:  `fwhms/fwhms_spec_cut.npz`, `regions-n.reg`
Output: `regions-n_[up,down].reg`

_Here, stop and check that the up/down regions look okay in ds9_
_Remove all dashes from filenames_

## Extract spectra (from multiple ObsIDs)

_Edit the shell scripts to use correct ObsIDs, weights, evt2 files!_

Code:   `../code-specs/[newregion,newspectrum,mergespectra].sh`
Input:  `regions_n_[up,down]_box.reg` (WCS fk5 coords)
Output: `spectra/[up,down]/*.[pi,rmf,arf]`

_Output filenames may require some manual adjustment, currently_

## Link spectra to backgrounds

_Run in debug mode first, check that filepaths to be linked are correct.
Then check in DS9 that the background/region linkings look correct.
Then, go ahead and modify files_

Code:   `../code/spec_linkbg.py`
Input:  spectra, `regions-n_[up,down].reg`, `backgrounds-n.reg`
Output: N/A (modifies spectra in place)

Also, yes, the background link is the ungrouped file.  See this
[CIAO thread](http://cxc.harvard.edu/ciao/threads/extended/index.html#up).

## Fit spectra to absorbed powerlaw with(out) Si line
Run with 32 bit python (`arch -i386 python`), `heainit` (not in same window as
CIAO), and run in same directory as spectra files (to resolve links)

Code:   `../code/spec_fit.py`
Input:  `spectra/[up,down]/*.[pi,rmf,arf]`
Output: `spectra/[up,down]/plot/plt_*.ps`,
        `spectra/[up,down]/fit/fit_*.[log,npz,json]`

Code:   `../code/spec_plot2pdf.sh`
Input:  `spectra/[up,down]/plot/plt_*.ps`
Output: `spectra/[up,down]/plot/plt_*.pdf`

NOTE: double check fit logs before proceeding, to make sure that background
spectra are linked correctly!

## Plot spectra and profiles

Code:   `../code-models/plotter_prfs_specs.ipynb`
Input:  `fwhms/fwhms.pkl`, `spectra/[up,down]/fit/fit_*.[log,npz,json]`
Output: `fwhms/plots/spec_*.pdf`

Uses both profiles and spectra to generate three-panel plot


Pipeline 3: fit FWHM data to filament models
============================================
Prerequisite: 1

## Add SNR information

Code:   `code-models/snr_catalog.py`

Add a method that initializes a `SupernovaRemnant` with appropriate constants.
Also stores a number of constants/magic numbers for fits


## Generate pre-cached table of FWHMs over parameter space
Let run overnight or so.  Clone generating script (basically your config file,
don't lose it) into `code-models` so it can import `models.py`.  When complete,
move outputs to `./tables/` and `chmod a-w [names]`.

For Tycho, it's set up to generate grids at fixed shock velocity (specify a few
different values and just run in separate terminal sessions)

Code:   `code-models/models.py`, called from `code-models/tables/gen-scripts/*`
Input:  bunch of magic numbers, settings, SNR object from `snr_catalog.py`
Output: `*_gen_2014-*_grid_*-*-*.[pkl,log,errlog]`

## Fit to simple/full/damped models
See Ressler et al. 2014 for exposition of models, calculations

Code:   `code-models/*-fit-tables.ipynb`
Input:  `./tables/*.pkl`, `data-*/regions-*/fwhms/fwhms_*.pkl`
Output: `data/fwhms/model-fits/*.[pkl, json]`

Make directories `model-fits` and `model-fits-damp` first.
Use IPython's parallel setup.  Generates interactive output in `.ipynb`,
absolutely necessary for debugging/etc.


Pipeline 4: generate tables for paper
=====================================
Prerequisite: 1, 2, 3, depends

## Data FWHMs, mE values, spectral fit parameters (all)

Code:   `../code-models/plotter_prfs_specs.ipynb`
Input:  `fwhms/fwhms.pkl`, `spectra/[up,down]/fit/fit_*.[log,npz,json]`
Output: table of FWHMs and mE values, table of absorbed power law fits,
        table of power law fits with gaussians for lines, or excised lines

Currently, fits accounting for lines aren't included in paper,
but the table is useful to get parameter values to cite in text




Example commands
================
Give a full pipeline thing when done (a single shell script would be nice).
For now this is a melange of text, not useful yet.

First generate data files and plots

    python ../../code/ds9projplotter.py -v \
        regions-n.reg profiles/prf \
        --pltroot profiles/plots/plt \
        --bands 0.7-1kev_mosaic.fits 1-1.7kev_mosaic.fits [...] \
        --labels '0.7-1kev' '1-1.7kev' [...]

Now generate data files for counts (no plots)

    python ../../code/ds9projplotter.py -v \
        regions-n.reg profiles/prf-cts \
        --bands 0.7-1kev_counts.fits 1-1.7kev_counts.fits [...] \
        --labels '0.7-1kev' '1-1.7kev' [...]

You can also generate, e.g. 3 band images, or subplot figures to better
compare various bands.  If data is already generated, you can omit --bands
argument so long as labels and profile root match filenames.

    python ../../code/ds9projplotter.py -v \
        regions-n.reg profiles/prf \
        --pltroot profiles/plots/plt-3band \
        --labels '0.7-1kev' '1-2kev' '2-7kev'

Subplots:

    python ../../code/ds9projplotter.py -v -s \
        regions-n.reg profiles/prf \
        --pltroot profiles/plots/plt-3band-sp \
        --labels '0.7-1kev' '1-2kev' '2-7kev'

.... break in continuity here

Now run a script to make spectra for them all, WITHOUT backgrounds.
Run the script 1x on the selected background regions too.

    # Use CIAO specextract commands
    ciao
    punlearn ardlib
    punlearn specextract
    pset specextract infile= \
    "../chandra/10095/repro/acisf10095_repro_evt2.fits[sky=@regions-good-3.ciaoreg]"
    pset specextract outroot="spectra/good-3/reg"
    pset specextract verbose=1
    specextract mode=h

    # Not implemented
    ciao
    python ciaoreg2spec.py -v 'regions.ciaoreg' 'spectra/reg'
    python ciaoreg2spec.py -v 'bkgs.ciaoreg' 'spectra/bkg'

    # Not implemented
    python update_spec_names_with_region_labels.py -v 'spectra/reg'

Link spectra to background regions. (must have CIAO initialized)

    ciao # if not initialized
    python spec_linkbg.py 'regions.ciaoreg' 'bkgs.ciaoreg' \
        'spectra/reg' 'spectra/bkg'



Example bash commands for pipeline
==================================
(out of date, at least as of regions-4-ext, or around early July 2014)

Some verbatim code for processing `regions-good-3`.
The relevant input and output files are:

    data/regions-good-3.reg                     # Hand-picked region file
    data/regions-good-3.ciaoreg
    data/backgrounds.reg                        # Hand-picked background file
    data/backgrounds.ciaoreg

    data/profiles/good-3/prf_{...}_{color}.dat  # Radial profile data
    plots/good-3/plt_{...}.png                  # Radial profile plots

    data/spectra/good-3/reg_{...}.{pi/arf/rmf}  # Spec + auxiliary
    data/spectra/good-3/reg_{...}.fitlog        # Fit parameters
    data/spectra/plots-good-3/plt_{...}.pdf     # Spec+model plots

Pipeline to produce this output (almost prepared to stitch together into a
single bash script like thing... but, need to fix a directory structure,
address CIAO/HEASOFT parameter file conflicts?, consider updates to scripts).
So for now, no need to chain together.  It's easy enough to run as-is.

    cd ~/Documents/snr-research/data

    # Generate regions
    ./ds9-tycho-rgb-unbin.sh

    # Plot and save radial profile data
    python ../code/ds9projplotter.py -v regions-good-3.reg \
        ../plots/good-3/plt -d profiles/good-3/prf

    # Convert to CIAO regions
    python ../code/ds9proj2box.py -v 2-7kev_mosaic.fits regions-good-3.reg
    regions-good-3.ciaoreg

    # Execute specextract
    ciao
    punlearn ardlib
    punlearn specextract
    pset specextract infile= \
        "../chandra/10095/repro/acisf10095_repro_evt2.fits \
        [sky=@regions-good-3.ciaoreg]"
    pset specextract outroot="spectra/good-3/reg"
    pset specextract verbose=1
    specextract mode=h

    # Link spectra with nearest background regions
    python ../code/spec_linkbg.py -v regions-good-3.ciaoreg backgrounds.ciaoreg
    spectra/good-3/reg spectra/bkg/bkg

    # Run XSPEC model fits, plot data and model, and save fit parameters
    heainit
    cd ./spectra/good-3/
    arch -i386 python ../../../code/spec_fitplot.py reg ../plots-good-3/plt
    # OR, with symlink...
    arch -i386 python spec_fitplot.py reg ../plots-good-3/plt

    # Convert ps to pdf
    spec_plot2pdf.sh  # Supply root "plt" when prompted

