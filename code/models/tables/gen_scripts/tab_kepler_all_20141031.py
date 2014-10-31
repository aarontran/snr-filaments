"""
Tabulate FWHMs for Kepler for
#1. no damping, all 6 mu values; default shock velocity (5000 km/s, 5 kpc)
#   20+ B0 values (ask for 20, expect ~30)
2. various magnetic damping scale lengths
   50+ B0 values (ask for 50, expect 10--70)

  `a_b` values I use are:
      0.5,    0.05,   0.04,   0.03,
      0.02,   0.015,  0.012,  0.011,
      0.01,   0.009,  0.008,  0.007,
      0.006,  0.005,  0.004,  0.003,  0.002

Default rminarc, Bmin = 5 microGauss.
151 eta2 values log-spaced between [0.01, 100] with 2x sampling in [0.1, 10].

Numbers now hardcoded into table generating script (makes more sense)
Call script with flag 0, 1, 2, 3 and it will load the correct values.
I don't know why I didn't do this before (since the principle is that this
script is essentially a "config" file of sorts).

to try to balance the load of calculations (I have gone diagonally in this
table).


Aaron Tran
2014 October 31
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np

import models
import snr_catalog as snrcat

AB_VALS_ALL = [[0.5,  0.015,  0.008, 0.003],
               [0.05, 0.012, 0.007, 0.006, 0.002],
               [0.04, 0.011, 0.01, 0.005],
               [0.03, 0.02,  0.009, 0.004]]

# Numbers are generated from regions-3, background subtracted
KEPLER_KEVS = np.array([0.7, 1.0, 2.0, 4.0])
KEPLER_DATA_MIN = np.array([ 2.0,  1.32, 1.20,  0.90])
KEPLER_DATA_MAX = np.array([ 3.3,  2.85, 2.65,  2.38])


def main():
    parser = argparse.ArgumentParser(description='Table generator')
    parser.add_argument('num', help=('Instance of the script to run; '
                                     'determines input a_b/v_s values. '
                                     'Value between 0-3'))

    args = parser.parse_args()
    idx = int(args.num)

    # vs/ab values to use in this instance
    abvals = AB_VALS_ALL[idx]

    # mu / eta2 / B0 grid settings
    mu_vals_damp = [1.]
    n_B0_damp = 50

    eta2_vals = [np.logspace(-2, -1, 25, base=10, endpoint=False),
                 np.logspace(-1, 1, 100, base=10, endpoint=False),
                 np.logspace(1, 2, 26, base=10, endpoint=True)]
    eta2_vals = np.sort(np.concatenate(eta2_vals))

    # Kepler settings (yes this is rather silly)
    snr = snrcat.make_kepler()

    kevs = KEPLER_KEVS
    data_min = KEPLER_DATA_MIN / 1.1
    data_max = KEPLER_DATA_MAX * 1.1

    # First make non-damped table
    if idx == 3:

        mu_vals_nodamp = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
        n_B0_nodamp = 20

        fname_nodamp = '{}_gen_{}_grid_{}-{}-{}.pkl'.format(snr.name,
                        datetime.now().strftime('%Y-%m-%d'), len(mu_vals_nodamp),
                        len(eta2_vals), n_B0_nodamp)

        models.maketab(snr, kevs, data_min, data_max,
                       mu_vals_nodamp, eta2_vals, n_B0_nodamp,
                       fname = fname_nodamp,
                       f_B0_init = 1.1, f_B0_step = 0.10,
                       idamp = False)

    # Then slew of damped tables
    for ab in abvals:

        fname_damp = '{}_gen_{}_grid_{}-{}-{}_ab-{:0.2e}.pkl'.format(snr.name,
                     datetime.now().strftime('%Y-%m-%d'),
                     len(mu_vals_damp), len(eta2_vals), n_B0_damp, ab)

        models.maketab(snr, kevs, data_min, data_max,
                       mu_vals_damp, eta2_vals, n_B0_damp,
                       fname = fname_damp,
                       f_B0_init = 1.1, f_B0_step = 0.10,
                       idamp = True, damp_ab = ab, damp_bmin = 5e-6)

if __name__ == '__main__':
    main()

