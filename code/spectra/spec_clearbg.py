"""
Script to clear BACKFILE entries from set of files
Aaron Tran
June 16, 2014

This is for when you run spec_linkbg.py and mess up the order of arguments
(and modify the wrong files...).  Modified from spec_linkbg.py

Initialize CIAO before running this script!
"""

import argparse
import os
import re

from ciao_contrib.runtool import dmhedit, dmhistory

def main():
    """Parse user input and update spectra background links
    KEY ASSUMPTION: numbering of region, background spectra
    MUST match numbering/ordering of region, background CIAOREG files
    """
    parser = argparse.ArgumentParser(description=
             'Clear BACKFILE entries for set of spectra')
    parser.add_argument('specroot', help='Directory stem for spectra')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    specroot = args.specroot
    verbose = args.verbose
    if verbose:
        print '\nFix spectra with stem: {}'.format(os.path.abspath(specroot))

    # Get number of spectra to update (count ungrouped spectrum files)
    pattern = os.path.basename(specroot) + r'_src[0-9]+\.pi'
    files = [f for f in os.listdir(os.path.dirname(specroot))
             if re.match(pattern, f)]
    print '{} files will be modified'.format(len(files)*2)

    # Set file header for each spectra (grouped and ungrouped)
    if verbose:
        print '\nUpdating file headers...'
    for num in xrange(len(files)):
        set_bkg(num+1, specroot, verbose) # Numbering starts at 1

    if verbose:
        print '\nDone!'


def set_bkg(num, rt, verbose=False):
    """Execute CIAO dmhedit to link spectra to backgrounds
    Check that specified root points to actual files
    (but I don't check that they are FITS files with spectra)
    
    Be careful -- this method modifies files
    (i.e., there is risk of data loss!)
    """
    # Create paths, with assumptions about filename structure
    reg_path = '{root}_src{number:d}.pi'.format(root=rt, number=num)
    reggrp_path = '{root}_src{number:d}_grp.pi'.format(root=rt, number=num)

    # Check for valid files
    if not (os.path.isfile(reg_path) and os.path.isfile(reggrp_path)):
        print 'One of these paths is bad:'
        print reg_path, reggrp_path
        raise Exception('ERROR: path does not exist!')
    if verbose:
        print '\nResetting BACKFILE for:'
        print ' {}\n {}'.format(reg_path, reggrp_path)

    # Set dmhedit parameters
    dmhedit.punlearn()
    dmhedit.filelist = 'none'
    dmhedit.operation = 'add'
    dmhedit.key = 'BACKFILE'
    dmhedit.value = 'none'
    if verbose:
        dmhedit.verbose = str(1)
    else:
        dmhedit.verbose = str(0)

    # Execute dmhedit on each spectrum
    for spec in [reg_path, reggrp_path]:
        dmhedit.infile = spec
        if verbose:
            print dmhedit()
        else:
            dmhedit()


if __name__ == '__main__':
    main()

