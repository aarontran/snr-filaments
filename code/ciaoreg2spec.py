"""
Script to generate spectra for CIAOREG region file
Aaron Tran
June 12, 2014

Currently: you MUST initialize CIAO before running this script!
"""

import argparse
from ciao_contrib.runtool import *

def main():
    """ Blah"""
    parser = argparse.ArgumentParser(description=
             'Generate individual spectra for regions in a CIAO region file.')
    parser.add_argument('ciaoreg', help='CIAO region file for specextract
                        stack')
    parser.add_argument('obsid_dir', help='Chandra ObsID directory (with repro/
                        from chandra_repro')
    parser.add_argument('outroot', help='Stem of output spectra')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    cr, obsdir, outroot = args.ciaoreg, args.obsid_dir, args.outroot
    verbose = args.verbose

    print 'Face it this is low priority for now.'
    print 'This is merely a wrapper for ciao\'s own functionality'
    print 'The next script down the line is important, however'

    ## Run the first time only, then don't reset
    # Set bad pixel file, evt2 file, default specextract parameters
    # TESTING WHETHER WE NEED THIS, IF WE SPECIFY
    # EVT FILE IN SPECEXTRACT
    # Nope.  Just punlearn ardlib first and you're good!!!
    punlearn ardlib
    punlearn specextract
    acis_set_ardlib acisf10095_repro_bpix1.fits verbose=1 absolutepath=yes
    pset specextract verbose=2

    # These are all set correctly by default, but just in case...
    pset specextract weight=yes
    pset specextract weightrmf=no
    pset specextract combine=no
    pset specextract bkgresp=yes

    $evt2 = "acisf10095_repro_evt2.fits"

    pset specextract infile="$evt2[sky=@profiles_good_cutback.ciaoreg]"
    pset specextract outroot="spectra/region-n"
    specextract mode=h  # Do not prompt


if __name__ == '__main__':
    main()
