"""
Tabulate FWHMs for SN1006

Idea is that whenever a new table is generated, a new config file should be
generated to record input parameters etc (and, how they were computed).

Aaron Tran
2014 September 2
"""

from __future__ import division

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
    # Set data
    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data_all = np.array([data for data, eps in SN1006_DATA.values()])
    data_min = np.amin(data_all/1.25, axis=0)
    data_max = np.amax(data_all*1.25, axis=0)

    # Set grid (values should be, preferably, sorted)
    mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
    eta2_vals = np.logspace(-2, 3, 60, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    # Set a few other twiddle-ables here
    # using default filename, SNR rminarc
    models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0,
                   f_B0_init = 1.1, f_B0_step = 0.10)


if __name__ == '__main__':
    main()
