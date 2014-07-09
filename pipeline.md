snr-research code pipeline
==========================
Aaron Tran
Summer 2014

N.B. this may not be fully up to date -- see `rsch-notes.md`, and talk to me
before using pipeline.

The hope is that it's straightforward enough to make sense out of
(see code documentation / command-line help).

The most fuzzy/hand-wavy part is region selection, see notes in data/ for that.
A lot of tweaking goes into 1. FWHM fitting, and 2. spectrum fitting.  So you
need to look at the code for that to see what numbers are being set / frozen /
thawed / whatever.

Below I give an abbreviated version of the pipeline, for quick reference and to
see where files go / what code does what.

Pipeline 1: radial profiles and fits
==================================

## Prerequisites
Generate energy band images, in both intensity flux units and in uncorrected
counts (no exposure/vignetting correction).  These will be used throughout.

## Region creation, cataloging
Open RGB image of the Tycho data.

    ./ds9-tycho-rgb.sh

Select regions of interest and backgrounds by hand.  Add text labels
(they will not show up in the CIAO regions).  Color code regions.
Pick regions with slight thermal uptick in back.  Save copies in fk5 and in
physical coordinates.

Output: `regions-n.reg`, `regions-n.physreg`

## Generate radial profiles from regions

Code:   `../code/ds9projplotter.py`
Input:  `regions-n.reg`
Output: `profiles/prf_[...].dat`, `profiles/prf-cts_[...].dat`

## Fit radial profiles and obtain FWHMs

Code:   `../code-profiles/profile_process.ipynb`
Input:  `regions-n.reg`, `regions-n.physreg`
Output: `fwhms/fwhm-fits.[txt, pkl, log]`

Pipeline 2: extract and fit spectra
=================================
Prerequisite: 1

## Prep/split regions based on profile fits

Code:   `../code/ds9projsplitter.py`
Input:  `fwhms/fwhm-fits.pkl`, `regions-n.reg`
Output: `regions-n-[up,down].reg`

Code:   `../code/ds9proj2ciao.py`
Input:  `regions-n-[up,down].reg`
Output: `regions-n-[up,down].ciaoreg`

## Extract spectra

Code:   CIAO specextract (Python wrapper `ciaoreg2spec.py` tbd)
Input:  `regions-n-[up,down].ciaoreg`
Output: `spectra/[up,down]/*.[pi,rmf,arf]`

*Repeat procedure to obtain background spectra, before continuing*

Code:   `../code/spec_linkbg.py`
Input:  `regions-n.ciaoreg`, `backgrounds-n.ciaoreg`
Output: N/A (modifies spectra in place)

## Fit spectra to absorbed powerlaw with Si line
Run with 32 bit python (`arch -i386 python`), `heainit` (not in same window as
CIAO), and run in same directory as spectra files (to resolve links)

Code:   `../code/spec_fit.py`
Input:  `spectra/[up,down]/*.[pi,rmf,arf]`
Output: `spectra/[up,down]/plots/plt_*.ps`,
        `spectra/[up,down]/fits/fit_*.[log,npz,json]`

Code:   `../code/spec_plot2pdf.sh`
Input:  `spectra/[up,down]/plots/plt_*.ps`
Output: `spectra/[up,down]/plots/plt_*.pdf`


Pipeline 3: fit FWHM data to filament models
============================================
Prerequisite: 1

## Fit to catastrophic dump transport model
Equation (6) of Ressler paper

Code:   `../code-models/fwhms_process.ipynb`

Still in development


Pipeline 4: generate plots for paper?
=====================================
Prerequisite: 1, 2 (so far)

## Plot spectra and profiles





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

