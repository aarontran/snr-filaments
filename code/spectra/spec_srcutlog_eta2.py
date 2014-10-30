"""
Parse srcutlog fit outputs from spec_fit.py, calculate eta2 as a function of mu

Print plaintext file of:
    srcut break frequencies (converted to linear space),
    eta2 values

Aaron Tran
Fall 2014
"""

from __future__ import division

import argparse
from datetime import datetime
from glob import glob
import json
import numpy as np
import os
import re

import json_utils as jsu

PLANCK_KEV_S = 4.135667517e-18  # Planck's constant; units: keV s
C_M = 1.82e18  # CGS from Sean's synchrotron relation \nu_m = c_m * E^2 * B
               # where C_M ~ 0.29 * 3*e/(4*pi*m^3*c^5)
REDCHI_WARN = 1.2

def main():
    """Parse input, run computations, spew output"""
    parser = argparse.ArgumentParser(description=
             ('Read in srcutlog fits from filestem; '
              'generate list of break frequencies (keV/h) and eta2 values'))

    parser.add_argument('fitroot', help='File stem of spectrum srcutlog fits')
    parser.add_argument('vsfile', help='Plaintext file of region shock velocities')
    parser.add_argument('outfile', help='Filename for output eta2 values')
    parser.add_argument('-m', '--mu-values', nargs='*', metavar=('mu'),
                        help=('mu values at which to calculate eta2 values. '
                              'Default: [0, 1/3, 0.5, 1, 1.5, 2]'))

    args = parser.parse_args()
    fitroot = args.fitroot
    vsfname = args.vsfile
    outfname = args.outfile
    mu_vals = args.mu_values
    if not mu_vals:  # [] or None
        mu_vals = [0, 1./3, 0.5, 1, 1.5, 2]
    else:
        mu_vals = map(float, mu_vals)

    vsvals = np.loadtxt(vsfname)

    # Require exact match to expected filename pattern, to be safe
    # ls_fits is only used to get number of files
    re_prog = re.compile(r'^' + fitroot + r'_src[0-9]+_grp.json$')
    ls_fits = glob(fitroot + '_src[0-9]*_grp.json')  # glob first
    ls_fits = [x for x in ls_fits if re_prog.match(x)]

    arr_out = []

    for n, vs in zip(xrange(1, len(ls_fits)+1), vsvals):  # vsvals ordered

        fitfname = fitroot + '_src{:d}_grp.json'.format(n)
        assert fitfname in ls_fits

        # Extract break frequency, check that chisqr is reasonable
        spec_info = jsu.load(fitfname)
        redchi = spec_info['fitstat'][1] / spec_info['dof']
        if redchi > REDCHI_WARN:  # Arbitrary
            _str = 'Warning: srcutlog fit for region {} has redchi {}'
            print _str.format(n,redchi)
        logbreak = spec_info['comps']['srcutlog']['break']['value']
        nu_break = 10**logbreak  # break freq in Hz
        en_break = PLANCK_KEV_S * nu_break  # break energy in keV

        row_out = [en_break]

        for mu in mu_vals:
            # Pull out all the (1+mu) powers and apply to nu_break/C_M first
            numer = (13.3)**4 * (2657)**(-(1-mu)) * (vs / 1e8)**4
            denom = (nu_break / (C_M * 100e-6))**(1+mu)
            eta2 = np.sqrt(numer/denom)
            row_out.append(eta2)

        arr_out.append(row_out)

    meta = 'srcutlog break frequencies and derived eta2 values\n'
    meta += 'Source spectrum fits: {}\n'.format(os.path.abspath(fitroot))
    meta += 'Shock velocity list file: {}\n'.format(os.path.abspath(vsfname))
    meta += 'eta2 values computed for mu values: {}\n'.format(mu_vals)
    meta += 'Generated: {}\n'.format(datetime.now())
    meta += 'Column 1 is break frequency (keV/h), rest are eta2 values'

    np.savetxt(outfname, arr_out, fmt='%.18f', header=meta)

    print 'Wrote srcutlog breaks and eta2 values to {}'.format(outfname)


if __name__ == '__main__':
    main()
