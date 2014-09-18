"""
Save intensity profiles and generate profile plots for a DS9 region file
containing only projection-type regions.
Generate profiles/plots for an arbitrary number of FITS files and labels.
Fig. 7 of Ressler et al. [ApJ, in press] illustrates what we're doing

Relies on pyds9 to extract 2D projection plots.
Uses ds9's RGB frames to force projection plots to pop-up and saves data from
the plot windows.

pyds9 documentation: http://hea-www.harvard.edu/RD/pyds9/
DS9 XPA reference: http://ds9.si.edu/ref/xpa.html

Aaron Tran
June 12, 2014
(last modified Sept 14, 2014)
"""

# TODO bug -- if data already exist, and code is asked to generate data/plots
# it throws an exception when it tries to generate plots.
# > Generating 22 plots...
# > First plot will hang a little longer
# > Traceback (most recent call last):
# >   File "../../code/ds9projplotter.py", line 186, in <module>
# >     main()
# >   File "../../code/ds9projplotter.py", line 87, in main
# >     generate_plots(rspecs, datroot, pltroot, labels, verbose, subplot)
# >   File "../../code/ds9projplotter.py", line 114, in generate_plots
# >     for flab, j in zip(labels, xrange(n)):
# > UnboundLocalError: local variable 'n' referenced before assignment

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

import ds9

from fplot import fplot, show_mplrc_settings
import regparse

def main():
    """Run from command line as, e.g.,
    python ds9projplotter.py -v '../data/profiles_all.reg' '../plots/plt'
    """
    parser = argparse.ArgumentParser(description=
             ('Plot and save data from DS9 projections. If only labels are '
              'given, it will just plot existing dat files'))

    parser.add_argument('infile', help='Input DS9 region file (fk5 coords)')
    parser.add_argument('datroot', help='Stem for output profile data')

    parser.add_argument('-p', '--pltroot', help='Stem for output plots')
    parser.add_argument('-f', '--files', help=('Files (.fits) to be processed,'
                                                ' one for each energy band'),
                        nargs='*')

    parser.add_argument('-l', '--labels', help='Band labels for plots/files',
                        nargs='+')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    parser.add_argument('-s', '--subplot', action='store_true',
                        help='Generate subplot plots')

    # Parse arguments
    args = parser.parse_args()
    regions_file = args.infile
    pltroot, datroot = args.pltroot, args.datroot
    bands, labels = args.files, args.labels  # NOTE bands = args.files...
    verbose, subplot = args.verbose, args.subplot

    # Check arguments for sanity
    # If -p is supplied, output plots
    # If -b is supplied, output data
    if bands is None:
        if not pltroot:
            parser.error('Nothing to be done!  Need bands or pltroot.')
        else:
            print 'No bands supplied, generating plots only'
    else:
        if len(bands) != len(labels):
            raise ValueError(('Non-zero number of bands '
                             'does not match number of labels!'))
        elif not pltroot:
            print 'No pltroot supplied, generating data only'

    # Parse projections
    rspecs = get_region_params(regions_file)
    if verbose:
        print 'Projections parsed from region file'

    # Generate data files, if desired
    if bands is not None:
        regparse.check_dir(datroot, verbose)  # only create if saving data
        generate_dat_files(rspecs, datroot, bands, labels)

    # Generate plots, if desired
    if pltroot:
        regparse.check_dir(pltroot, verbose)
        generate_plots(rspecs, datroot, pltroot, labels, verbose, subplot)

    if verbose:
        print 'Done!'


########################
# Various helper methods
########################

def generate_plots(rspecs, datroot, pltroot, labels,
                   verbose=False, subplot=False):
    """Iterate through generated dat files and make multi-band plots"""

    if subplot:
        n = len(labels)
        plt.figure(figsize=(4*n,4))
    else:
        plt.figure(figsize=(6,4))

    if verbose:
        print 'Generating {} plots...'.format(len(rspecs))
        print 'First plot will hang a little longer'

    # Make a plot for each region
    for i in xrange(len(rspecs)):
        # Plot line/points for each label
        for flab, j in zip(labels, xrange(n)):
            dat_fname = '{0}_{1:02d}_band_{2}.dat'.format(datroot, i+1, flab)
            a = np.loadtxt(dat_fname)
            if subplot:
                plt.subplot(1, n, j+1)  # j indexes labels
                plt.plot(a[:,0], a[:,1], '-o')  # Use default color cycle
                plt.title('Region {:g}, band {}'.format(i+1, flab))
                fplot('Radial dist. (?)', 'Intensity? (?)')
            else:
                plt.plot(a[:,0], a[:,1], '-o')  # Use default color cycle
        # Format overall plot
        if not subplot:
            fplot('Radial dist. (?)', 'Intensity? (?)')
            plt.title('Region %g' % (i+1))
            plt.legend(labels, loc='best')
        plt.tight_layout()
        plt.savefig('{0}_{1:02d}.png'.format(pltroot, i+1), dpi=150)
        plt.clf()
        if verbose:
            print 'Saved: {0}_{1:02d}.png'.format(pltroot, i+1)


def generate_dat_files(rspecs, datroot, bands, labels):
    """Start DS9 and iterate through bands/regions to save data"""
    d = ds9.ds9()
    d.set('rgb')
    d.set('rgb red')

    # Save plaintext projection data
    # Idea: minimize file (band) loading operations
    for fname, flab in zip(bands, labels):
        d.set('file ' + fname)  # Load a band
        for i in xrange(len(rspecs)):
            d.set('regions', rspecs[i])  # Load a region
            d.set('rgb red')  # Plot projection data
            dat_fname = '{0}_{1:02d}_band_{2}.dat'.format(datroot, i+1, flab)
            d.set('plot {0} save {1}'.format(d.get('plot'), dat_fname))
            d.set('regions delete all')
    d.set('exit')


def get_region_params(fname):
    """Read in DS9 region file (version 4.1), return list of XPA arguments
    for projections ONLY!
    """
    # Read strings from the region file
    with open(fname, 'r') as f:
        lines = f.readlines()
        wcs = lines[2][:-1]  # Coordinate system TODO bad variable name
        if lines[0] != '# Region file format: DS9 version 4.1\n':
            print 'Warning: potentially invalid region file!'
            print 'First line was: ' + lines[0]
        if wcs != 'fk5':
            raise Exception('Regions must be in sky (fk5) coordinates; got ' +
                            wcs + 'instead')

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
