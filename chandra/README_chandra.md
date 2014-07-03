README for raw Chandra observations
===================================

Note: ObsID folders (and other stuff in this folder) are not tracked by git.
Chandra data is downloaded from [WebChaSeR](http://cda.harvard.edu/chaser/).

ObsID folders `10[0-9]+` contain raw and reprocessed (via `chandra_repro`)
data from Hughes' 750 ks observation of Tycho's supernova remnant.
Observations were made ~April 2009 using ACIS-I.
The exposure times are as follows:

* 10093, 118.35 ks
* 10094, 89.97 ks
* 10095, 173.37 ks
* 10096, 105.72 ks
* 10097, 107.43 ks
* 10902, 39.53 ks
* 10903, 23.92 ks
* 10904, 34.7 ks
* 10906, 41.12 ks

ObsID folders `671*` and `7366` contain raw data from Reynolds' 750 ks
observation of Kepler's supernova remnant.
Observations were made ~April to August 2006 using ACIS-S (I presume S3 chip).
The exposure times are as follows:

* 6714, 157.82 ks (repro 2014 July 2)
* 6715, 159.13 ks (repro 2014 July 2)
* 6716, 158.02 ks
* 6717, 108.0 ks
* 6718, 110.0 ks
* 7366, 52.0 ks (repro 2014 July 1)

Reprocessing 6714, 6715 gave the following error (3, 4 times respectively):

    # acis_process_events (CIAO 4.6): The following error occurred 3 times:
    dsAPEPULSEHEIGHTERR -- WARNING: pulse height is less than split threshold when performing serial CTI adjustment.

Attempt to reprocess 6716 errored out immediately with:

    # chandra_repro (20 Nov 2013): ERROR CRC check failed 0x184fcdfe != 0x43ff83b4L 
