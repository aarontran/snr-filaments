"""
Script to split ds9 regions into inner/outer pieces (split into two files)

Aaron Tran
2014 July 8

This script is intended for our SNR radial profiles (Tycho, Kepler, etc.).

To compute spectra, select inner and outer regions delineated by
1. rim between largest measured FWHMs.  Find lowest downstream and highest
   upstream FWHM locations (yes, may mix limits from different bands).
2. downstream FWHM to most downstream local minimum, i.e. boundary of fit domain
   for profile fitting scripts.  Again, checks cut location for all bands
   and selects the most conservative cut (largest possible region).

One input region file will be split into two for specextract/XSPEC chain.
Each file contains only inner or only outer pieces, to preserve numbering.

Script makes MANY assumptions about structure of region dictionary.
So, if dictionary structure is changed... update accordingly.
"""

import argparse
import cPickle as pickle
import numpy as np
import os
import re

import regparse


def main():
    """Generate split regions from cmd line arguments"""
    parser = argparse.ArgumentParser(description=
             ('Split DS9 region file of projections into inner/outer pieces, '
              'according to profile fit cuts and FWHMs'))
    parser.add_argument('fitpkl', help=('Pickled region dictionary from '
                        'profile fitting scripts'))
    parser.add_argument('fitreg', help='DS9 fk5 region file of projections')
    parser.add_argument('regimg', help='Valid FITS file for input regions')
    parser.add_argument('-o', '--outreg', help=('Root name for output region '
                        'files.  Inferred from input region file by default'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose mode')

    # Parse arguments
    args = parser.parse_args()
    f_pkl, f_reg, f_img = args.fitpkl, args.fitreg, args.regimg
    f_oreg = args.outreg
    verbose = args.verbose

    if f_oreg is None:
        f_oreg = os.path.splitext(f_reg)[0]
    if verbose:
        print 'Output root is {}'.format(f_oreg)

    with open(f_pkl, 'r') as f:
        regions = pickle.load(f)

    # Convert input region file to physical coordinates
    f_reg_phys = f_reg + 'tmp'
    if verbose:
        print ('Converting input to physical coordinates; '
               'tmp file {}'.format(f_reg_phys))
    regparse.conv_fk5_to_phys(f_reg, f_img, f_reg_phys)
    rspecs, headers = regparse.load_ds9reg(f_reg_phys, aux=True)
    os.remove(f_reg_phys)  # Clean-up

    rspecs_down = []
    rspecs_up = []

    # Go through each region, extract correct cuts
    # and save new subsetted regions

    # NOTE Ordering/numbering of regions must match!  This is NOT ENFORCED.
    if verbose:
        print 'Applying FWHM/fit cuts from {}'.format(f_pkl)
    for rspec, n in zip(rspecs, regions.keys()):
        reg = regions[n]
        x1, x2, x3 = scale_cuts(reg['info']['spec_cuts'], reg)

        rspecs_down.append(subset_proj(rspec, x1, x2))
        rspecs_up.append(subset_proj(rspec, x2, x3))

    # Save as ds9, physical coords, then resave as wcs,fk5 coords
    f_down = '{}_down.reg'.format(f_oreg)
    f_up = '{}_up.reg'.format(f_oreg)
    if verbose:
        print 'Writing split regions to temporary physical coord files'
    regparse.write_ds9reg(f_down+'tmp', rspecs_down, headers)
    regparse.write_ds9reg(f_up+'tmp', rspecs_up, headers)

    # Convert to ds9, fk5 coordinates
    regparse.conv_reg_coords(f_down+'tmp', f_img, f_down,
                             fmt='ds9', sys='wcs', sky='fk5')
    regparse.conv_reg_coords(f_up+'tmp', f_img, f_up,
                             fmt='ds9', sys='wcs', sky='fk5')
    # Cleanup
    os.remove(f_down+'tmp')
    os.remove(f_up+'tmp')

    if verbose:
        print 'Wrote output: {}, {}'.format(f_down, f_up)
        print 'Done!'


def scale_cuts(x, reg):
    """Convert cuts from arcsec, measured from first data point (x=0),
    to fractional position along projection length"""
    px2as = reg['info']['px2arcsec']
    length = reg['info']['length']  # Pixels, not floored
    return x/(px2as*length)


def subset_proj(rstr, a, b):
    """Generates projection subsetted along radial distance between a, b
    where a, b \in [0,1] and a < b
    Input, output: ds9 region strings in physical coordinates
    """
    rtype, rnums, rprops = regparse.regparse(rstr)
    x1, y1, x2, y2, t = rnums

    x1a = x1 + a*(x2-x1)
    y1a = y1 + a*(y2-y1)

    x2b = x1 + b*(x2-x1)
    y2b = y1 + b*(y2-y1)

    return '{0}({1},{2},{3},{4},{5}){6}'.format(rtype,x1a,y1a,x2b,y2b,t,rprops)


if __name__ == '__main__':
    main()

