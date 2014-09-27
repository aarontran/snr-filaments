"""
Script to link region spectra with nearest backgrounds
Aaron Tran
2014 July 12
(last modified: 2014 Sept 26)

Initialize CIAO before running this script!
"""

import argparse
import os
import re
import sys

from ciao_contrib.runtool import dmhedit, dmhistory

# NOTE we must tell CIAO python where pyds9 is located
sys.path.append('/usr/local/lib/python2.7/site-packages')
import regparse

def main():
    """Parse user input and update spectra background links
    KEY ASSUMPTION: numbering of region, background spectra
    MUST match numbering/ordering of region, background region files
    """
    parser = argparse.ArgumentParser(description=
             ('Link spectra to closest background spectra (modifies FITS '
              'header keywords). Acts on orig + grouped files. '
              'Directory stem for spectra is whatever precedes _src{:d}*.*'))

    parser.add_argument('regfile', help=('DS9 fk5 region file for spectra to '
                                         'be linked (regions must be boxes!)'))
    parser.add_argument('bkgfile', help='DS9 fk5 region file for backgrounds')
    parser.add_argument('imgfile', help='FITS file w/ image for regions')
    parser.add_argument('regroot', help='Directory stem for region spectra')
    parser.add_argument('bkgroot', help='Directory stem for bkg spectra')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='debug mode (no files modified)')

    args = parser.parse_args()
    regfile, bkgfile = args.regfile, args.bkgfile
    imgfile = args.imgfile
    regroot, bkgroot = args.regroot, args.bkgroot
    verbose, debug = args.verbose, args.debug

    if verbose:
        print 'Regions and backgrounds from files:'
        print ' Regions: {}'.format(os.path.abspath(regfile))
        print ' Backgrounds: {}'.format(os.path.abspath(bkgfile))
        print '\nLinking spectra with filename roots:'
        print ' Regions: {}'.format(os.path.abspath(regroot))
        print ' Backgrounds: {}'.format(os.path.abspath(bkgroot))

    # Convert ds9 fk5 to CIAO physical coords -- temporary files
    if verbose:
        print '\nConverting DS9 fk5 to CIAO physical coords'
    regfile_c = regfile + '.tmp.ciaoreg'
    bkgfile_c = bkgfile + '.tmp.ciaoreg'
    regparse.conv_reg_coords(regfile, imgfile, regfile_c,
                             fmt='ciao', sys='physical')
    regparse.conv_reg_coords(bkgfile, imgfile, bkgfile_c,
                             fmt='ciao', sys='physical')

    # CIAO numbering starts at 1, not 0, hence the key numbering...
    # Centers, used as keys, should be unique.  But, if not unique,
    # then doesn't mtter which number is used
    if verbose:
        print '\nParsing center coordinates'
    reg_ctrs = get_centers(regfile_c, verbose)
    bkg_ctrs = get_centers(bkgfile_c, verbose)
    reg_dict = {reg_ctrs[x]:(x+1) for x in xrange(len(reg_ctrs))}
    bkg_dict = {bkg_ctrs[x]:(x+1) for x in xrange(len(bkg_ctrs))}

    # Remove CIAO files -- only needed them to get region strings/centers
    if verbose:
        print '\nCleaning up temporary CIAO region files'
    os.remove(regfile_c)
    os.remove(bkgfile_c)

    # Find closest bkg to each region
    if verbose:
        print '\nMatching regions to nearest backgrounds...'
    rb_dict = {}
    for ctr in reg_dict.keys():
        bkg_ctr = closest(ctr, bkg_dict.keys()) # Closest background center
        rb_dict[reg_dict[ctr]] = bkg_dict[bkg_ctr]  # Keys, values are numbers

    # Assign bkgfile to each region (spec + grouped spec)
    if verbose:
        print '\nUpdating file headers...'
    # Yes, bad practice, but this lets us load/use pyds9 before loading CIAO
    #sp = subprocess.Popen(['/bin/bash', '-i', '-c', 'ciao'])#call('ciao', shell=True)
    #sp.communicate()
    #from ciao_contrib.runtool import dmhedit, dmhistory
    for r in rb_dict.keys():
        set_bkg(r, rb_dict[r], regroot, bkgroot, verbose, debug)

    if verbose:
        print '\nDone!'


def set_bkg(rnum, bnum, r_rt, b_rt, verbose=False, debug=False):
    """Execute CIAO dmhedit to link spectra to backgrounds
    Check that r_rt and b_rt point to actual files
    (but I don't check that they are actual FITS files with spectra)

    Be careful -- this method modifies files
    (i.e., there is risk of data loss!)
    """
    # Create paths, with assumptions about filename structure
    reg_path = '{root}_src{number:d}.pi'.format(root=r_rt, number=rnum)
    reggrp_path = '{root}_src{number:d}_grp.pi'.format(root=r_rt, number=rnum)
    bkg_path = '{root}_src{number:d}.pi'.format(root=b_rt, number=bnum)

    # Check for valid files
    if not (os.path.isfile(reg_path) and os.path.isfile(reggrp_path)
            and os.path.isfile(bkg_path)):
        print 'One of these paths is bad:'
        print reg_path, reggrp_path, bkg_path
        raise Exception('ERROR: path does not exist!')

    # Now, find RELATIVE path from reg/reggrp files, to bkg file
    reg2bkg_path = os.path.relpath(bkg_path, os.path.dirname(reg_path))
    if verbose:
        print '\nSetting file headers for:'
        print ' {}\n {}'.format(reg_path, reggrp_path)
        print ' Relative path to bkg: {}'.format(reg2bkg_path)

    # Set dmhedit parameters
    dmhedit.punlearn()
    dmhedit.filelist = 'none'
    dmhedit.operation = 'add'
    dmhedit.key = 'BACKFILE'
    dmhedit.value = '\'{}\''.format(reg2bkg_path) # Single quotes for paths
    if verbose:
        dmhedit.verbose = str(1)
    else:
        dmhedit.verbose = str(0)

    # Execute dmhedit on each spectrum
    for spec in [reg_path, reggrp_path]:
        dmhedit.infile = spec
        if not debug:
            if verbose:
                print dmhedit()
            else:
                dmhedit()
        else:
            print 'DEBUG: A call to dmhedit would occur here'


def get_centers(fname, verbose=False):
    """Given a CIAO region file, parse and return sensible region centers
    Regions must be rectangles, boxes, circles, or ellipses.
    REGION ARITHMETIC IS NOT SUPPORTED -- ONLY ISOLATED REGIONS.

    Input: fname (str) is a conforming CIAO region filename
    Output: list of region centers (2-tuple of floats, physical coords)
    """
    with open(fname, 'r') as f:
        # Check region file type
        header = f.readline()
        if 'CIAO version 1.0' not in header:
            print 'Bad header: {}'.format(header)
            raise Exception('ERROR: invalid region file!')
        # Parse each region
        centers = []
        for line in f:
            centers.append(get_reg_center(line))

    return centers


def get_reg_center(regstr):
    """Parse a single region string, obtain center for valid regions
    Regions must be rectangles, boxes, circles, or ellipses.
    REGION ARITHMETIC IS NOT SUPPORTED -- ONLY ISOLATED REGIONS.
    Future updates may deal with more complex regions

    Input: regstr (str) is a single string from a CIAO region file
    Output: 2-tuple of floats, physical coords of region's center
    """
    # List of regions with spec: name(xcenter,ycenter,...)
    # and, where (x,y)center actually correspond to a geometric centroid
    # of a 2-D region (so annulus, sector, pie, point are not okay)
    # http://cxc.harvard.edu/ciao/ahelp/dmregions.html
    reglist = ['rect', 'rectangle', 'box', 'rotbox', 'cir', 'circle',
               'ell', 'ellipse']

    specs = re.split('\s*[\(\)\,]\s*', regstr)
    reg = specs[0]
    if reg.lower() not in reglist:
        print 'Incompatible region: {}'.format(reg)
        raise Exception('ERROR: invalid region for center parsing')

    return (float(specs[1]), float(specs[2]))


########################
# Point geometry methods
########################

def closest(a, pts):
    """Find the closest point in pts to point a"""
    c = pts[0]
    sqdist_ac = sqdist(a, c)
    for p in pts[1:]:
        sqdist_ap = sqdist(a, p)
        if sqdist_ap < sqdist_ac:
            c = p
            sqdist_ac = sqdist_ap
    return c


def sqdist(a, b):
    """Square of Euclidean distance between two 2-tuples of floats
    Using square to avoid a square root computation"""
    x1, y1 = a
    x2, y2 = b
    return (x1-x2)**2 + (y1-y2)**2


if __name__ == '__main__':
    main()

