README for raw Chandra observations
===================================

Tycho's SNR observation
-----------------------

Note: ObsID folders (and other stuff in this folder) are not tracked by git.
Chandra data is downloaded from [WebChaSeR](http://cda.harvard.edu/chaser/).

ObsID folders `10[0-9]+` contain raw and reprocessed (via `chandra_repro`)
data from Hughes' 750 ks observation of Tycho's supernova remnant.
Observations were made ~April 2009 using ACIS-I.
The exposure times are as follows:

* 10093, 118.35 ks
* 10094,  89.97 ks
* 10095, 173.37 ks
* 10096, 105.72 ks
* 10097, 107.43 ks
* 10902,  39.53 ks
* 10903,  23.92 ks
* 10904,  34.7  ks
* 10906,  41.12 ks

(total exposure time: 734.11 ks, if copied/summed correctly)

The folder `imgs_tycho` contains `output[0-9]+\.fits` files, which are
relatively small (~5 MB) image files for each observation.  Obtained from Nina
Coyle, who obtained them from Brian Williams -- used by scripts extracting
spectra from all ObsIDs and merging them together.

To try reproducing Brian's image mosaics, use xygrid="3300:4900:1,3300:4900:1"
(see procedure outlined below for Kepler).

Folder `tycho_thermal` just holds CIAO `merge_obs` output for 1.7-2.5 keV
image, for poster session.

Kepler's SNR observation
------------------------

ObsID folders `671*` and `7366` contain raw data from Reynolds' 750 ks
observation of Kepler's supernova remnant.
Observations were made ~April to August 2006 using ACIS-S (I presume S3 chip).
The exposure times are as follows:

* 6714, 157.82 ks (repro 2014 July 2)
* 6715, 159.13 ks (repro 2014 July 2)
* 6716, 158.02 ks (repro 2014 July 11)
* 6717, 108.0  ks (repro 2014 July 11)
* 6718, 110.0  ks (repro 2014 July 11)
* 7366,  52.0  ks (repro 2014 July 1)

Reprojected images in `out_kepler/` generated with CIAO 4.6 as:

    ciao
    punlearn ardlib
    punlearn reproject_obs
    reproject_obs infiles="6714,6715,6716,6717,6718,7366" outroot="kepler_reproj/"

> Warning: the merged event file kepler\_reproj/merged\_evt.fits
>    should not be used to create ARF/RMF/exposure maps because
>       the RA_NOM keyword varies by 0.00990496859004 (limit is 0.0003)
>       the DEC_NOM keyword varies by 0.066592311011 (limit is 0.0003)
>       the ROLL_NOM keyword varies by 176.679053066 (limit is 1.0)

Flux images in `kepler/` generated with CIAO 4.6 (manually pick xygrid) as:

    punlearn flux_obs
    pset flux_obs bands="0.7:1:0.85, 1:2:1.5, 2:7:4.5"
    pset flux_obs xygrid="3640:4540:1,3630:4530:1"
    flux_obs kepler_reproj/ kepler/

Count images (not created) can be generated as:

    punlearn dmcopy
    dmcopy "merged_evt.fits[EVENTS][energy=1000:1700][bin x=3640:4540:1,y=3630:4530:1]" \
           outfile="1-1.7kev_counts.fits"

varying the energy range and outfile appropriately.  This procedure for
generating count images follows an email from Brian J. Williams.


Cassiopeia A observation
------------------------

ObsID folders 4634-4639, 5196, 5319, 5320 contain raw data from Una Hwang's 1 Ms
observation of Cassiopeia A.
Observations were made February to May 2004 using ACIS-S
The exposure times are as follows:

* 4634, 148.62 ks (repro 2014 July 12)
* 4635, 135.04 ks (repro 2014 July 12)
* 4636, 143.48 ks (repro 2014 July 14)
* 4637, 163.5  ks (repro 2014 July 14)
* 4638, 164.53 ks (repro 2014 July 14)
* 4639,  79.05 ks (repro 2014 July 14)
* 5196,  49.52 ks (repro 2014 July 14)
* 5319,  42.26 ks (repro 2014 July 14)
* 5320,  54.37 ks (repro 2014 July 14)


