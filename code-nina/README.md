Notes on Nina/Brian's shell scripts for merging multiple CIAO specextract spectra
-----------
(2014 July 10)

File directory setup Nina uses is as follows:

    data/   NAMEcoord.reg
            newregion.sh
            insample.lis
            outsample.lis
            
            spectra/    newspectrum.sh      # Nina's version (works on Mac)
                        newspectrum_or.sh   # Brian's older version
                        
                        aspfile.lis         # *.lis files only needed for Brian's version
                        bpixfile.lis        # Nina's version calls specextract for each ObsID
                        mskfile.lis         # so no need
                        pbkfile.lis

newregion.sh
------------
Open ds9 window (ds9&) + call heainit before running
Edit paths to CHANDRA repro files here + in insample/outsample.lis
Call as:

./newregion.sh NAME

newspectrum.sh
--------------
Call CIAO before running newspectrum.sh
edit the filenames to CHANDRA repro files here as well
Call as:

./newspectrum.sh NAME