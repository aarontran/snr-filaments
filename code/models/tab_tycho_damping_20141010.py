"""
Tabulate FWHMs for Tycho for various magnetic damping scale lengths
Uses default shock velocity (3.6e8 * 3/2.3 = 4.696e8 cm/s ...)
Uses default damping Bmin (5 microGauss)
Using default rminarc

151 eta2 values log-spaced between [0.01, 100] with 2x sampling in [0.1, 10].
~50 B0 values (ask for 30), sampled to give ~linear spacing in FWHM widths

data_min/data_max set to a bit smaller range than before

`a_b` values I use are:

    0.5,    0.05,   0.04,   0.03,
    0.02,   0.01,   0.009,  0.008,
    0.007,  0.006,  0.005,  0.004

But, feed damping numbers into scripts as:

    python tab_tycho_damping_141010.py  0.5   0.01   0.005
    python tab_tycho_damping_141010.py  0.05  0.009  0.004
    python tab_tycho_damping_141010.py  0.04  0.008  0.007
    python tab_tycho_damping_141010.py  0.03  0.02   0.006

to try to balance the load of calculations (I have gone diagonally in this
table).


Aaron Tran
2014 October 10
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np

import models
import snr_catalog as snrcat

# Numbers are generated from regions-4 with simple 2-exponential fit
# Max observed FWHM for 0.7-1keV was 6.092
# But I use 10 to be more consistent w/ general trend
# (matches previous tables)
TYCHO_KEVS = np.array([0.7, 1.0, 2.0, 3.0, 4.5])
TYCHO_DATA_MIN = np.array([1.628, 1.685, 1.510, 1.465, 1.370])
TYCHO_DATA_MAX = np.array([10.0, 8.866, 6.901, 7.508, 5.763])

def main():
    # Read damping length
    parser = argparse.ArgumentParser(description='Table generator')
    parser.add_argument('abvals', help='Damping lengths (%% of shock radius)',
                        nargs='*')
    args = parser.parse_args()
    abvals = map(float, args.abvals)

    # Set SNR parameters
    snr = snrcat.make_tycho()

    # Set data
    kevs = TYCHO_KEVS
    data_min = TYCHO_DATA_MIN / 1.1
    data_max = TYCHO_DATA_MAX * 1.1

    # Set grid (values should be, preferably, sorted)
    mu_vals = [1.]
    eta2_vals = [np.logspace(-2, -1, 25, base=10, endpoint=False),
                 np.logspace(-1, 1, 100, base=10, endpoint=False),
                 np.logspace(1, 2, 26, base=10, endpoint=True)]
    eta2_vals = np.sort(np.concatenate(eta2_vals))
    n_B0 = 30  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing (esp. for small n_B0)

    for ab in abvals:

        # Set output filename base (end with .pkl)
        # Using custom name with damping length ab appended, for Tycho
        fname = '{}_gen_{}_grid_{}-{}-{}_ab-{:0.2e}.pkl'.format(snr.name,
                datetime.now().strftime('%Y-%m-%d'),
                len(mu_vals), len(eta2_vals), n_B0, ab)

        # Set a few other twiddle-ables here
        models.maketab(snr, kevs, data_min, data_max,
                       mu_vals, eta2_vals, n_B0,
                       fname = fname,
                       f_B0_init = 1.1, f_B0_step = 0.10,
                       idamp = True, damp_ab = ab)

if __name__ == '__main__':
    main()

