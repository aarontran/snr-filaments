"""
Tabulate FWHMs for Tycho, given a value for shock velocity

Idea is that whenever a new table is generated, a new config file should be
generated to record input parameters etc (and, how they were computed).

Inefficient and ugly... but it should do the job?

PART 2 -- resuming where I hit box-length errors due to rminarc too small for
highest energy band, when diffusion is large.
Two separate calls to models.maketab (for each mu) to account for slightly
different eta2 gridding, attempting to save time.

Aaron Tran
2014 August 1
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np

import snr_catalog as snrcat
import models_all as models

# Numbers are generated from regions-4 with simple 2-exponential fit
# Max observed FWHM for 0.7-1keV is actually 6.092
# But I use 10 to be more consistent w/ general trend, and to be more
# conservative in tabulation (esp. since some regions don't have FWHM at
# 0.7-1kev)
TYCHO_KEVS = np.array([0.7, 1.0, 2.0, 3.0, 4.5])
TYCHO_DATA_MIN = np.array([1.628, 1.685, 1.510, 1.465, 1.370])
TYCHO_DATA_MAX = np.array([10.0, 8.866, 6.901, 7.508, 5.763])


def main():
    # Read shock velocity
    parser = argparse.ArgumentParser(description='Table generator')
    parser.add_argument('vs', help='Shock velocity (cm/s)')
    args = parser.parse_args()
    vs = float(args.vs)

    # Set SNR parameters
    snr = snrcat.make_tycho()
    snr.vs = vs

    # Set data
    kevs = TYCHO_KEVS
    data_min = TYCHO_DATA_MIN / 1.25
    data_max = TYCHO_DATA_MAX * 1.25

    # Set default grid
    eta2_vals_orig = np.logspace(-2, 2, 50, base=10)
    eta2_vals_orig = np.sort(np.append(eta2_vals_orig, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    # Set rminarc for gridding (determined manually) (NEW LARGER VALUES)
    rminarc = [18.5, 15., 10.77, 10.77, 10.77]
    # approx ~ TYCHO_DATA_MAX * 1.25 * f_rminarc, f_rminarc = 1.25

    # ------------------------------------
    mu_vals = [1.5]
    eta2_vals = eta2_vals_orig[eta2_vals_orig > 15]
    fname = '{}_gen_2014-07-30_grid_6-100-20_vs-{:0.2e}_part2_mu-1.50.pkl'.format(
            snr.name, vs)

    models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0,
                   rminarc = rminarc,
                   fname = fname,
                   f_rminarc = float('NaN'),  # Irrelevant, just being explicit
                   f_B0_init = 1.1,
                   f_B0_step = 0.10)

    # ------------------------------------
    mu_vals = [2]
    eta2_vals = eta2_vals_orig
    fname = '{}_gen_2014-07-30_grid_6-100-20_vs-{:0.2e}_part2_mu-2.00.pkl'.format(
            snr.name, vs)

    models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0,
                   rminarc = rminarc,
                   fname = fname,
                   f_rminarc = float('NaN'),  # Irrelevant, just being explicit
                   f_B0_init = 1.1,
                   f_B0_step = 0.10)


if __name__ == '__main__':
    main()
