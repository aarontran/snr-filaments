"""
Tabulate FWHMs for SN1006 for various magnetic damping scale lengths
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

    python ....py  0.5   0.01   0.005
    python ....py  0.05  0.009  0.004
    python ....py  0.04  0.008  0.007
    python ....py  0.03  0.02   0.006

to try to balance the load of calculations (I have gone diagonally in this
table).


Aaron Tran
2014 October 13
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np

import models
import snr_catalog as snrcat

# SN 1006 numbers (from Sean's script, originally)
SN1006_KEVS = np.array([0.7, 1.0, 2.0])
SN1006_DATA = {}
SN1006_DATA[1] = np.array([35.5, 31.94, 25.34]), np.array([1.73, .97, 1.71])
SN1006_DATA[2] = np.array([23.02, 17.46, 15.3]), np.array([.35,.139, .559])
SN1006_DATA[3] = np.array([49.14, 42.76,29.34]), np.array([1.5, .718, .767])
SN1006_DATA[4] = np.array([29, 23.9, 16.6]), np.array([.9, .39, .45])
SN1006_DATA[5] = np.array([33.75, 27.2, 24.75 ]), np.array([2.37,.62,.61])

def main():
    # Read damping length
    parser = argparse.ArgumentParser(description='Table generator')
    parser.add_argument('abvals', help='Damping lengths (%% of shock radius)',
                        nargs='*')
    args = parser.parse_args()
    abvals = map(float, args.abvals)

    # Set SNR parameters
    snr = snrcat.make_SN1006()

    # Set data
    kevs = SN1006_KEVS
    data_all = np.array([data for data, eps in SN1006_DATA.values()])
    data_min = np.amin(data_all/1.1, axis=0)
    data_max = np.amax(data_all*1.1, axis=0)

    # Set grid (values should be, preferably, sorted)
    mu_vals = [1.]
    eta2_vals = [np.logspace(-2, -1, 25, base=10, endpoint=False),
                 np.logspace(-1, 1, 100, base=10, endpoint=False),
                 np.logspace(1, 2, 26, base=10, endpoint=True)]
    eta2_vals = np.sort(np.concatenate(eta2_vals))
    n_B0 = 50  # In practice, you'll usually get ~1.5 to 2x as many points
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


