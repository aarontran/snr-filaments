README for data (mosaics, regions, etc)
=======================================
Aaron Tran
2014 June 11

(last modified: 2014 June 26, but still rather out of date)

Data/images
-----------
Brian's mosaics of the ~750 ks Hughes observation, 2009.

    0.7-1kev_mosaic.fits
    0.7-1kev_mosaic_bin.fits
    1-2kev_mosaic.fits
    1-2kev_mosaic_bin.fits
    2-7kev_mosaic.fits
    2-7kev_mosaic_bin.fits
    2009_3to7kev.fits

Additional energy bands (breaking 2-7 keV apart)

    2-3kev_mosaic.fits
    2-4kev_mosaic.fits
    3-4p5kev_mosaic.fits
    4-7kev_mosaic.fits
    4p5-7kev_mosaic.fits

Single, 3 GB file with all event counts

    merged_evt.fits

Count files for all mosaics (no exposure/vignetting correction), generated from
`merged_evt.fits`.

    *_counts.fits

Raw observations, downloaded from the CXC Webchaser interface.

    chandra/

Data analysis stuff
-------------------

Notes on region selection

    profiles_notes.md

Short script to load ds9 with the mosaicked, binned data for region selection.

    ds9-tycho-rgb.sh

All regions, good, bad, ugly.

    profiles_all.ciaoreg
    profiles_all.reg
    profiles_all_pic.jpeg

Subset of all regions, just the good ones

    profiles_good.ciaoreg
    profiles_good.reg
    profiles_good_bg.ciaoreg
    profiles_good_bg.reg
    profiles_good_cutback.ciaoreg
    profiles_good_cutback.reg
    profiles_good_pic.jpeg

Regions for testing

    test-regs


