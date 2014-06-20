"""
Plot intensity profiles from a three band image (you pick the bands)
by parsing a ds9 region file with projection-type regions.
Also saves profile data if so desired

Relies on pyds9 to extract 2D projection plots.
This will probably be replaced with a better program/approach, since
1. ds9's projection is a bit of a blackbox
2. code is rather slow since it has to plot data, write/read data every time
   (I couldn't find a way to just get numbers/values from plots directly)
but, it does work.

Fig. 7 of Ressler et al. [ApJ, in review/press?] illustrates what we're doing

pyds9 documentation: http://hea-www.harvard.edu/RD/pyds9/
DS9 XPA reference: http://ds9.si.edu/ref/xpa.html

Aaron Tran
June 12, 2014
(last modified June 17, 2014)
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

from fplot import fplot, show_mplrc_settings
import ds9


def main():
    """Run from command line as, e.g.,
    python ds9projplotter.py -v '../data/profiles_all.reg' '../plots/plt'
    """
    parser = argparse.ArgumentParser(description=
             'Plot profiles from DS9 region files with projections')

    parser.add_argument('infile', help='Input DS9 region file (fk5 coords)')
    parser.add_argument('outroot', help='Stem of output plots')
    # RGB files for profiles
    parser.add_argument('-r', '--red', help='red band FITS file',
        default='../data/0.7-1kev_mosaic_unbin.fits')
    parser.add_argument('-g', '--green', help='green band FITS file',
        default='../data/1-2kev_mosaic_unbin.fits')
    parser.add_argument('-b', '--blue', help='blue band FITS file',
        default='../data/2-7kev_mosaic_unbin.fits')
    # Labels for legends
    parser.add_argument('--rlabel', help='red band label',
        default='0.7-1 keV')
    parser.add_argument('--glabel', help='green band label',
        default='1-2 keV')
    parser.add_argument('--blabel', help='blue band label',
        default='2-7 keV')
    # Additional flags
    parser.add_argument('-d', '--data',
                        help='Save projection data to directory stem')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    
    # Parse arguments
    args = parser.parse_args()
    rfits, gfits, bfits = args.red, args.green, args.blue
    regions_file = args.infile
    outroot = args.outroot
    datroot = args.data
    verbose = args.verbose

    # Check that outroot, datroot are valid paths
    outdir = os.path.dirname(outroot)
    if outdir is not '' and not os.path.isdir(outdir):
        print 'Outroot {} in nonexistent directory'.format(outroot)
        print 'Creating directory {}'.format(outdir)
        os.makedirs(outdir)
    if datroot:
        datdir = os.path.dirname(datroot)
        if not os.path.isdir(datdir):
            print 'Data root {} in nonexistent directory'.format(datroot)
            print 'Creating directory {}'.format(datdir)
            os.makedirs(datdir)

    # Start DS9 with RGB image
    d = ds9.ds9()
    d.set('rgb')
    d.set('rgb red')
    d.set('file ' + rfits)
    d.set('rgb green')
    d.set('file ' + gfits)
    d.set('rgb blue')
    d.set('file ' + bfits)
    d.set('rgb lock scale yes')
    d.set('scale asinh')

    # Parse projections
    rspecs = get_region_params(regions_file)
    if verbose:
        print 'Projections parsed from region file'

    # Generate plots
    if verbose:
        print 'Generating {} plots...'.format(len(rspecs))
        print 'First plot will hang a little longer'
    plt.figure(figsize=(6,4))

    # Load regions one-by-one in DS9
    for i in xrange(len(rspecs)):
        d.set('regions', rspecs[i])

        # Cycle through RGB colors to get different plot data
        for clr in ['red', 'green', 'blue']:
            d.set('rgb ' + clr)
            # Save the data somewhere if so desired
            if datroot:
                dat_fname = '{0}_{1:02d}_{2}.dat'.format(datroot, i+1, clr)
            else:
                dat_fname = '{}_temp.dat'.format(clr)
            d.set('plot {0} save {1}'.format(d.get('plot'), dat_fname))
            a = np.loadtxt(dat_fname)
            plt.plot(a[:,0], a[:,1], '-o'+clr[0])

        # Now that all bands are plotted, format and save plot
        fplot('Radial dist. (?)', 'Intensity? (?)')
        plt.title('Region %g' % (i+1))
        plt.legend(('0.7-1 keV', '1-2 keV', '2-7 keV'), loc='best')
        #plt.tight_layout()
        plt.savefig('{0}_{1:02d}.png'.format(outroot, i+1), dpi=150)
        plt.clf()
        d.set('regions delete all')

    # Clean up
    if not datroot:
        if verbose:
            print 'Deleting temporary files'
        for clr in ['red', 'green', 'blue']:
            os.remove(clr + '_temp.dat')

    d.set('exit')
    if verbose:
        print 'Done!'


def get_region_params(fname):
    """Read in DS9 region file (version 4.1), return list of XPA arguments
    for projections ONLY!
    """
    # Read strings from the region file
    with open(fname, 'r') as f:
        lines = f.readlines()
        wcs = lines[2][:-1]  # Coordinate system
        if lines[0] != '# Region file format: DS9 version 4.1\n':
            print 'Warning: potentially invalid region file!'
            print 'First line was: ' + lines[0]
        if wcs != 'fk5':
            print 'Warning: potentially problematic coordinate system!'
            print 'May work fine, just not tested'
            print 'Coord. syst. is: ' + lines[2][:-1]

    # Manipulate each string and save to list
    lines = filter(lambda x: '# projection' in x, lines[3:])
    projspecs = []
    for ln in lines:
        lnsplit = ln[2:-1].split(')')  # Remove leading '#', trailing '\n'
        if lnsplit[1] == '':
            r = lnsplit[0] + ')'  # No optional arguments
        else:
            r = lnsplit[0] + ') #' + lnsplit[1] # Add octothorpe
        projspecs.append('%s; %s' % (wcs, r))

    return projspecs


if __name__ == '__main__':
    main()
