README for data (mosaics, regions, etc)
=======================================
Aaron Tran
2014 June 11

NB very out of date -- mainly documents files/images stored here

Data/images
-----------
Brian's mosaics of the ~750 ks Hughes observation, 2009.

    0.7-1kev_mosaic.fits
    1-2kev_mosaic.fits
    2-7kev_mosaic.fits

Additional energy bands (breaking 2-7 keV apart)

    2-3kev_mosaic.fits
    2-4kev_mosaic.fits
    3-4p5kev_mosaic.fits
    4-7kev_mosaic.fits
    4.5-7kev_mosaic.fits

Single, 3 GB file with all event counts

    merged_evt.fits

Count files for all mosaics (no exposure/vignetting correction), generated from
`merged_evt.fits` as:

    dmcopy "merged_evt.fits[EVENTS][energy=700:1000][bin x=3300:4900:1,y=3300:4900:1]" 0.7-1kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=1000:2000][bin x=3300:4900:1,y=3300:4900:1]" 1-2kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=1000:1700][bin x=3300:4900:1,y=3300:4900:1]" 1-1.7kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=2000:7000][bin x=3300:4900:1,y=3300:4900:1]" 2-7kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=2000:3000][bin x=3300:4900:1,y=3300:4900:1]" 2-3kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=2000:4000][bin x=3300:4900:1,y=3300:4900:1]" 2-4kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=3000:4500][bin x=3300:4900:1,y=3300:4900:1]" 3-4.5kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=4000:7000][bin x=3300:4900:1,y=3300:4900:1]" 4-7kev_counts.fits
    dmcopy "merged_evt.fits[EVENTS][energy=4500:7000][bin x=3300:4900:1,y=3300:4900:1]" 4.5-7kev_counts.fits

`TYCHO_IF1.fits` -- VLA radio image, 1994, from Stephen Reynolds (in turn, from
David Moffett?).  Presumably project AM0437.  From header: units Jy/beam, freq
1.375 GHz, A config.  Not sure if 1.375, 1.625 GHz freq observations could have
been somehow combined or anything.

`ds9-tycho-rgb-*.sh` -- convenience scripts to load RGB images


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


