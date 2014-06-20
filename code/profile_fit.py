"""
Script to attempt profile fits
See Ressler et al. [in review, ApJ]

Aaron Tran
June 17, 2014
(last modified June 17, 2014)
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

from fplot import fplot, show_mplrc_settings


def main():
    """Run from command line as, e.g.,
    python ds9projplotter.py -v '../data/profiles_all.reg' '../plots/plt'
    """
    parser = argparse.ArgumentParser(description=
             'Plot profiles from DS9 region files with projections')

    parser.add_argument('infile', help='Input DS9 region file (fk5 coords)')
    parser.add_argument('outroot', help='Stem of output plots')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    
    # Parse arguments
    args = parser.parse_args()
    verbose = args.verbose

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

    if verbose:
        print 'Done!'


if __name__ == '__main__':
    main()
