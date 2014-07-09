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

Region creation, cataloging
---------------------------
Open an RGB image of the Tycho data.

    ./ds9-tycho-rgb.sh

Select the regions of interest and backgrounds by hand.  Add text labels
(they will not show up in the CIAO regions).  Color code regions.

    data/regions-all.reg
    data/regions-n/regions-n.reg
    data/bkg-n/bkg-n.reg

In general -- pick regions either (1) with a bit of thermal uptick in the back,
for fitting, or (2) stringently avoiding said thermal emission.  You need both,
at the end.


Generate radial profiles from regions
-------------------------------------

Make plots of radial profiles and save profiles to plaintext.  Do this twice --
once for profiles in intensity units (with exposure/vignetting corrections),
once for profiles in count units.  E.g., as follows

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


Fit radial profiles and obtain FWHMs
------------------------------------
Use iPython notebook (at shell: `ipython notebook`)
`code-profiles/profile_process.ipynb`.  Supply a few configuration arguments
(e.g., labels to fit/store together, same labels given to `ds9projplotter.py`)
and notebook will output various files (pkl, txt) containing fit information.

This information is needed to then investigate FWHM-energy dependence (and
hence B field amplification / diffusion models), as well as generate XSPEC
spectra to check for thermal contamination (fit domains and FWHMs delineate
where to obtain spectra).


Region manipulation (spectra)
-----------------------------
Commands executed from data/regions-n/

Convert ds9 region files to CIAO region files
(for both regions of interest and backgrounds)

    python ../../code/ds9proj2ciao.py -v ../2-7kev_mosaic.fits \
        regions-n.reg regions-n.ciaoreg

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

Plot spectra with fits to absorbed powerlaw (`phabs*powerlaw`), and output
fit parameters to same folder as spectra (can change later).
1. Must run with 32 bit Python
2. Don't run in same window as CIAO for now (conflict not addressed)
3. must run in same directory as spectra files

    heainit
    arch -i386 python ../code/spec_fitplot.py \
                      spectra/good-2/reg spectra/plots-good-2/plt

File organization:

    data/spectra/bkg/bkg_src*[_grp].[pi,arf,rmf]
    data/spectra/good-2/reg_ ...
    data/spectra/test-2/reg_ ...



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

