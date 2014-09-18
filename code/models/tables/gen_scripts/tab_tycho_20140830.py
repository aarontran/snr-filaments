"""
Tabulate FWHMs for Tycho, given a value for shock velocity

Idea is that whenever a new table is generated, a new config file should be
generated to record input parameters etc (and, how they were computed).

Previously used shock velocities:
4.21e8, 4.52e8, 4.83e8, 5.14e8

Min, max shock velocities from regions-4 assuming d = 3 kpc:
4.50e8, 5.20e8

Today, let's use:
[4.59e8, 4.76e8, 4.94e8, 5.11e8]
which gives an even spacing around our currently assumed velocities.

Aaron Tran
2014 August 30
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np

import models
import snr_catalog as snrcat

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

    # Set grid (values should be, preferably, sorted)
    mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
    eta2_vals = np.logspace(-2, 3, 60, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    # Set rminarc for gridding (determined manually)
    rminarc = np.array([18.5, 15., 10.77, 10.77, 10.77]) # from 20140730 tycho pt. 2
    # approx ~ TYCHO_DATA_MAX * 1.25 * f_rminarc, f_rminarc = 1.25

    # Set output filename base (end with .pkl)
    # Using custom name with shock velocity appended, for Tycho
    fname = '{}_gen_{}_grid_{}-{}-{}_vs-{:0.2e}.pkl'.format(snr.name,
            datetime.now().strftime('%Y-%m-%d'), len(mu_vals), len(eta2_vals),
            n_B0, vs)

    # Set a few other twiddle-ables here
    models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0,
                   fname = fname,
                   f_B0_init = 1.1,
                   f_B0_step = 0.10,
                   rminarc = rminarc)


if __name__ == '__main__':
    main()
