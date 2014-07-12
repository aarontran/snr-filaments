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

Reprocessing gives this error sometimes (with variable # of occurrences)

    # acis_process_events (CIAO 4.6): The following error occurred 3 times:
    dsAPEPULSEHEIGHTERR -- WARNING: pulse height is less than split threshold when performing serial CTI adjustment.

Resolution -- no need to worry.  From CIAO ahelp on `acis_process_events`:
> When the CTI adjustment is applied to events on the back-illuminated CCDs
> (ACIS-S1 and S3), sometimes one of the pulse heights in a 3x3 pixel event
> island can drop below the split threshold if it was above the threshold
> before the adjustment. In the end, pixels that are below the split threshold
> are ignored when the total pulse height and energy are computed.
>
> If the number of times is small, then the warning may be safely ignored.

From `merge_obs`:

> Warning: the merged event file `out_kepler/merged_evt.fits`
> should not be used to create ARF/RMF/exposure maps because
>   the RA\_NOM keyword varies by 0.00990496859004 (limit is 0.0003)
>   the DEC\_NOM keyword varies by 0.066592311011 (limit is 0.0003)
>   the ROLL\_NOM keyword varies by 176.679053066 (limit is 1.0)


Cassiopeia A observation
------------------------

ObsID folders 4634-4639, 5196, 5319, 5320 contain raw data from Una Hwang's 1 Ms
observation of Cassiopeia A.
Observations were made February to May 2004 using ACIS-S
The exposure times are as follows:

* 4634, 148.62 ks (repro 2014 July 12)
* 4635, 135.04 ks (repro 2014 July 12)
* 4636, 143.48 ks
* 4637, 163.5  ks
* 4638, 164.53 ks
* 4639,  79.05 ks
* 5196,  49.52 ks
* 5319,  42.26 ks
* 5320,  54.37 ks


