README, Kepler data (fits, regions, etc)
========================================

Image files:

    0.7-1kev_counts.fits
    0.7-1kev_mosaic.fits
    1-2kev_counts.fits
    1-2kev_mosaic.fits
    2-7kev_counts.fits
    2-7kev_mosaic.fits

Generated as usual (see `snr-research/chandra/`).
Counts files generated from `merged_evts.fits` as:

    dmcopy "merged_evt.fits[EVENTS][energy=1000:2000][bin x=3640:4540:1,y=3630:4530:1]" outfile="1-2kev_counts.fits"
    dmcopy "merged_evt.fits[EVENTS][energy=700:1000][bin x=3640:4540:1,y=3630:4530:1]" outfile="2-7kev_counts.fits"
    dmcopy "merged_evt.fits[EVENTS][energy=2000:7000][bin x=3640:4540:1,y=3630:4530:1]" outfile="2-7kev_counts.fits"

`ds9-kepler-rgb.sh` -- shell script to load region files into DS9 with
convenient scale parameters

Region selections:

* `regions-1/`: first round, relatively narrow regions around Kepler SE limb
* `bkg-1/`: background circles drawn for `regions-1/`, all equally sized
