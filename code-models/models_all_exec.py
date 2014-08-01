"""
Code to execute model fits for measured FWHMs (using models, tables, etc
given in models_all.py

Here is where we should generate tables, figures, check effects of various
parameters, et cetera.  The other code should do the heavy lifting and
implementation, without worrying about details of SNRs.

This is the code where I should generate plots of chi-squared space, check that
the tables I generate are well sampled, call table generating functions, et
cetera.

Aaron Tran
2014 July 26
"""

from __future__ import division

import cPickle as pickle
from datetime import datetime
import lmfit
import matplotlib.pyplot as plt
import numpy as np
from numpy import f2py
import scipy as sp
from scipy import optimize
import sys

from fplot import fplot
import snr_catalog as snrcat
import models_all as models

# Could throw this into another short script? remember to call
# fm.readfglists() before anything else though.
import fullmodel as fm
fm.readfglists()  # Only do this once

# TEMPORARY (FOR TESTING PURPOSES ONLY)
SN1006_KEVS = np.array([0.7, 1.0, 2.0])
SN1006_DATA = {}
SN1006_DATA[1] = np.array([35.5, 31.94, 25.34]), np.array([1.73, .97, 1.71])
SN1006_DATA[2] = np.array([23.02, 17.46, 15.3]), np.array([.35,.139, .559])
SN1006_DATA[3] = np.array([49.14, 42.76,29.34]), np.array([1.5, .718, .767])
SN1006_DATA[4] = np.array([29, 23.9, 16.6]), np.array([.9, .39, .45])
SN1006_DATA[5] = np.array([33.75, 27.2, 24.75 ]), np.array([2.37,.62,.61])

# TEMPORARY (FOR TESTING PURPOSES, FIRST TABULATION ONLY)
TYCHO_KEVS = np.array([0.7, 1.0, 2.0, 3.0, 4.5])
TYCHO_DATA_MIN = np.array([1.628, 1.685, 1.510, 1.465, 1.370])
TYCHO_DATA_MAX = np.array([10.0, 8.866, 6.901, 7.508, 5.763])

# Note -- the max observed for 0.7-1keV is actually 6.092, for regions-4
# BUT, to be consistent with the general trend seen in mins/maxes, and to be
# more conservative in tabulation, I arbitrarily choose a larger value of 10

# These numbers are generated from regions-4, with simple 2-exponential fit
# using the following short script:
#import cPickle as pickle
#import numpy as np
#with open('fwhms.pkl', 'r') as fpkl:
#    data = pickle.load(fpkl)
#bands = ['1-1.7kev', '3-4.5kev', '2-3kev', '0.7-1kev', '4.5-7kev']
#for b in bands:
#    fwhm_min = np.nanmin([data[n][b]['meas']['fwhm'] for n in data])
#    fwhm_max = np.nanmax([data[n][b]['meas']['fwhm'] for n in data])
#    print 'Band {} min/max are {}, {}'.format(b, fwhm_min, fwhm_max)


def main():
    """main stuff.  Edit this at will to generate desired output..."""


def vistable(tab_fname):
    """Make plot of grid, using tab output from tab_eta2().
    Only shows grid for one (arbitrary, first) value of mu
    """
    with open(tab_fname, 'r') as fpkl:
        tab = pickle.load(fpkl)

    mu_vals = tab.keys()
    mu_dict = tab[mu_vals[0]]
    for eta2 in mu_dict.keys():
        B0_vals = mu_dict[eta2][0]  # [1] for FWHM values!
        plt.plot(eta2*np.ones(len(B0_vals)), B0_vals, 'ob')

    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$\eta_2$'); plt.ylabel(r'$B_0$')
    plt.title(r'\eta_2, B_0 grid for \mu = {}'.format(mu_vals[0]))
    plt.show()


def generate_SN1006_table(mu):
    """What it says.  Configured for single mu value.
    Move to a configuration script if any more SN1006 tables are generated."""

    #mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
    #mu_vals = np.array([1./3, 1./2, 1])  # Kolmogorov, Kraichnan, Bohm
    mu_vals = [mu]
    eta2_vals = np.logspace(-2, 2, 50, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data/1.51, axis=0)
    data_max = np.amax(data*1.51, axis=0)

    fname = 'sn1006_grid-1-100-20_2014-07-28_mu-{:0.2f}_partest.pkl'.format(mu)

    tab = models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals,
                         n_B0, fname=fname)

    return tab


# ======================================
# Some ad hoc functions to check effects
# of varying certain parameters, or just
# spewing out some numbers (for now...)
# ======================================


def make_table7_SN1006():
    """Make numbers blagh
    This reproduces Table 7 of Ressler et al.
    """

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS

    for flmt in [1,2,3,4,5]:
        print 'Filament {}'.format(flmt)
        print 'mu\teta2\t\t\tB0'
        for mu in [0, 1/3, 1/2, 1, 1.5, 2]:
            data, eps = SN1006_DATA[flmt]
            res = simple_fit(snr, kevs, data, eps, mu)

            print '{:0.2f}\t{:0.2f} +/- {:0.2f}\t\t{:0.2f} +/- {:0.2f}'.format(
                    res.values['mu'],
                    res.values['eta2'], res.params['eta2'].stderr,
                    res.values['B0']*1e6, res.params['B0'].stderr*1e6)

            # Next, error from confidence intervals...
            # Error bars on B0 and eta2 will be tremendous in many cases.


def check_effects():
    """Make numbers blagh
    What happens when we change
    compression ratio (i.e. v0) ?

    And, what happens when we change shock velocity (vs and v0 together)?
    I think the result will be, B0 changes dramatically...
    Given what I see with compression ratio.
    """
    mu = 1
    compratios = [4, 6, 8, 20]

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS

    for cr in compratios:

        snr.v0 = snr.vs/cr
        print '\ncompression ratio = {}'.format(cr)

        for flmt in [1,2,3,4,5]:

            data, eps = SN1006_DATA[flmt]
            res = simple_fit(snr, kevs, data, eps, mu)

            print ('Filament {}: mu = {:0.2f}\teta2 = {:0.2f} +/- {:0.2f}\t'
                   'B0 = {:0.2f} +/- {:0.2f}').format(flmt, res.values['mu'],
                    res.values['eta2'], res.params['eta2'].stderr,
                    res.values['B0']*1e6, res.params['B0'].stderr*1e6)

def check_tycho_effects():
    """A TEMPORARY METHOD... better version later?

    Checking effects of shock velocity on B0 and eta2 values
    Also, may as well check effect of distance
    """

    mu = 1
    shock_velocs = np.array([1.92e8, 3.11e8, 3.58e8, 4.06e8]) * 3 / 2.3
    # From Williams et al. 2013, scaled for D = 3 kpc
    # Values are min, min of {vs|vs>2500}, mean of {vs|vs>2500}, max of {...}
    # (I am trying to ignore the regions of anomalously low shock velocity)

    snr = snrcat.make_tycho()
    kevs = TYCHO_KEVS
    # Region 1, from regions-4 set on Tycho, simple 2-exp fit
    data = np.array([2.72, 2.2, 2.22, 1.88, 1.90])
    eps = np.array([0.23, 0.05, 0.08, 0.08, 0.11])

    np.set_printoptions(precision=2)
    print 'Simple model fits, Tycho, region 1'
    print 'Using mu = {}, data: {}'.format(mu, data)
    print 'Investigating effect of changing shock velocity'

    for vs in shock_velocs:
        # Currently, I am computing tycho radius w D=3kpc assumption
        # So we are only looking at vs variation, ignoring dist
        snr.vs = vs
        snr.v0 = snr.vs/snr.cratio
        res = models.simple_fit(snr, kevs, data, eps, mu)

        print ('Shock velocity = {:0.2e}: mu = {:0.2f}\t'
               'eta2 = {:0.2f} +/- {:0.2f}\t'
               'B0 = {:0.2f} +/- {:0.2f}').format(vs, res.values['mu'],
                res.values['eta2'], res.params['eta2'].stderr,
                res.values['B0']*1e6, res.params['B0'].stderr*1e6)


    print '\nNow investigating remnant distance effect'
    dists = np.array([2.3, 3, 4, 5])

    for d in dists:
        snr.vs = 3.58e8 * d / 2.3
        snr.v0 = snr.vs / snr.cratio
        snr.rs = 1.077e19 * d / 3
        res = models.simple_fit(snr, kevs, data, eps, mu)

        print ('Dremnant = {} (vs = 3.58e8 * d/2.3): mu = {:0.2f}\t'
               'eta2 = {:0.2f} +/- {:0.2f}\t'
               'B0 = {:0.2f} +/- {:0.2f}').format(d, res.values['mu'],
                res.values['eta2'], res.params['eta2'].stderr,
                res.values['B0']*1e6, res.params['B0'].stderr*1e6)


if __name__ == '__main__':
    main()
