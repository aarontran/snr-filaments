snr-research code pipeline
==========================
Aaron Tran
Summer 2014


Region creation, cataloging
---------------------------

Open an RGB image of the Tycho data (binned or unbinned).

    ./ds9-tycho-rgb.sh
    ./ds9-tycho-rgb-unbin.sh

Select the regions of interest and backgrounds by hand.  Add text labels
(they will not show up in the CIAO regions).  Color code regions.
File organization:

    data/regions-all.reg
    data/regions-good-2.reg
    data/backgrounds.reg


Region manipulation (spectra)
-----------------------------

Convert ds9 region files to CIAO region files
(for both regions of interest and backgrounds)

    python ds9proj2box.py -v \
        ../data/2-7kev_mosaic.fits \
        ../data/profiles_good_cutback.reg \
        ../data/profiles_good_cutback.ciaoreg

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

    # Not implemented yet
    ciao
    python ciaoreg2spec.py -v 'regions.ciaoreg' 'spectra/reg'
    python ciaoreg2spec.py -v 'bkgs.ciaoreg' 'spectra/bkg'

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

Region manipulation (radial profiles)
-------------------------------------

Make plots of radial profiles and save profiles to plaintext

    python ds9projplotter.py -v \
        ../data/regions-good-2.reg \
        ../plots/good-2/plt \
        -d ../data/profiles/good-2/prf

Here's the newer version now optimized (not really) for an arbitrary number of
energy bands... first generate data files for everything, then output n-band
plots using the same script.

    python ../code/ds9projplotter.py -v regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-all -b 0.7-1kev_mosaic.fits 1-2kev_mosaic.fits 2-3kev_mosaic.fits 2-4kev_mosaic.fits 2-7kev_mosaic.fits 3-4p5kev_mosaic.fits 4-7kev_mosaic.fits 4p5-7kev_mosaic.fits -l '0.7-1kev' '1-2kev' '2-3kev' '2-4kev' '2-7kev' '3-4p5kev' '4-7kev' '4p5-7kev'
    
	# For overlaid plots
    python ../code/ds9projplotter.py -v regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-3band -l '0.7-1kev' '1-2kev' '2-7kev'
    python ../code/ds9projplotter.py -v regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-4band -l '0.7-1kev' '1-2kev' '2-4kev' '4-7kev'
    python ../code/ds9projplotter.py -v regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-5band -l '0.7-1kev' '1-2kev' '2-3kev' '3-4p5kev' '4p5-7kev'
	
	# For subplots
	python ../code/ds9projplotter.py -sv regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-3band-sp -l '0.7-1kev' '1-2kev' '2-7kev'
	python ../code/ds9projplotter.py -sv regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-4band-sp -l '0.7-1kev' '1-2kev' '2-4kev' '4-7kev'
	python ../code/ds9projplotter.py -sv regions-good-3-allback.reg profiles/good-3-allback/prf -p profiles/good-3-allback/plt-5band-sp -l '0.7-1kev' '1-2kev' '2-3kev' '3-4p5kev' '4p5-7kev'

Geez.  To get 3/4/5 band subplot images for only one region, run:
`open plt-*band-sp_02.png` (replacing 02 with number of region of choice).
And it will make it super easy to compare regions in multiple bands.


Apply fit model to spectra using `profile_fits.ipynb` (under development).


Example bash commands for pipeline
==================================

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

