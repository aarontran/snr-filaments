Notes on Nina/Brian's shell scripts for merging multiple CIAO specextract spectra
-----------
(2014 July 10)

File directory setup for scripts, out of the box, is:

    data/   [NAME]coord.reg
            newregion.sh
            insample.lis
            outsample.lis
            
            spectra/    newspectrum.sh      # Nina's version (works on Mac)
                        newspectrum_or.sh   # Brian's older version
                        
                        aspfile.lis   # *.lis files for Brian's older version
                        bpixfile.lis  # Nina's version calls specextract for
                        mskfile.lis   # each ObsID; no need for list files
                        pbkfile.lis

In all files (`newregion.sh`, `newspectrum.sh`, `*sample.lis`), edit paths to
Chandra repro/ files as needed.

newregion.sh
------------
Open ds9 window (ds9&) + call heainit before running script as:

    ./newregion.sh [NAME]

The original region file should be in CIAO WCS (fk5) coordinates;
this script generates output CIAO physical coordinates for EACH ObsID.
Due to slightly different quirks of each observation, these coordinates are in
general NOT THE SAME -- hence this script.

newspectrum.sh
--------------
Call heainit and CIAO before running script as:

    ./newspectrum.sh [NAME]

Binning
-------
Lastly, after generating combined spectra with .wrmf, .warf files, use HEASOFT
utility `grppha` to bin spectra if needed.  From email from Brian to Nina:

> ... it will prompt you for an input file (your total "spectrum.pi" file) and
> an output name (give it a name like "spectrum\_grp.pi"). At the next line,
> simply type
>
> `group min 15`
>
> then type exit.
