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
    parser.add_argument('regimg', help='Valid FITS file for input region file')
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

    # Manipulate DS9 region specs in PHYSICAL coordinates
    f_reg_phys = f_reg + '.tmpphys'
    regparse.conv_fk5_to_phys(f_reg, f_img, f_reg_phys)
    rspecs, headers = regparse.load_ds9reg(f_reg_phys, aux=True)
    os.remove(f_reg_phys)

    rspecs_down = []
    rspecs_up = []
    
    # Numbering should match order of region file
    for rspec, n in zip(rspecs, regions.keys()):
        reg = regions[n]
        x1, x2, x3 = get_cuts(reg)
        x1, x2, x3 = scale_cuts([x1,x2,x3], reg)

        rspecs_down.append(subset_proj(rspec, x1, x2))
        rspecs_up.append(subset_proj(rspec, x2, x3))

    # Save to new files
    regparse.write_ds9reg('{}-down.physreg'.format(f_oreg), rspecs_down, headers)
    regparse.write_ds9reg('{}-up.physreg'.format(f_oreg), rspecs_up, headers)

    if verbose:
        print 'Done!'


def scale_cuts(x, reg):
    """Convert cuts from arcsec, measured from first data point (x=0),
    to fractional position along projection length"""
    px2as = reg['info']['px2arcsec']
    length = reg['info']['length']  # Pixels, not floored
    return x/(px2as*length)


def get_cuts(reg):
    """Cut locations for spectra in arcsec, parsed from region dictionary"""
    labels = reg.keys()
    labels.remove('info')

    # Problem: what about bad FWHMs? -- need to manually blacklist
    x_btw = np.nanmin([reg[lab]['meas']['fwhm-lims'][0] for lab in labels])
    x_max = np.nanmax([reg[lab]['meas']['fwhm-lims'][1] for lab in labels])

    x_min = -1
    for lab in labels:
        # Only use cuts where FWHM could be fitted
        if np.isfinite(reg[lab]['meas']['fwhm']):
            ind = reg[lab]['cut']
            x_cut = reg[lab]['data'][0][ind]
            if x_min > x_cut or x_min == -1:
                x_min = x_cut
    
    return x_min, x_btw, x_max


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

