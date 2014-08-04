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

#import cPickle as pickle
#from datetime import datetime
import lmfit
import matplotlib.pyplot as plt
import numpy as np
#from numpy import f2py
#import scipy as sp
#from scipy import optimize
#import sys

from fplot import fplot
#import fullmodel as fm
#fm.readfglists()
#import snr_catalog as snrcat
import models_all as models



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


# Table from Press et al., Section 14.5 (1st ed.),
# 15.6 (3rd ed.); see also Avni, Y. (1976, ApJ)
# Would be good to compute properly in a function
CI_2_CHISQR = {}
CI_2_CHISQR[0.683] = 1.0
CI_2_CHISQR[0.90] = 2.71
CI_2_CHISQR[0.9545] = 4.0
CI_2_CHISQR[0.99] = 6.63


def main():
    """main stuff.  Edit this at will to generate desired output..."""


# =====================================================
# Parameter space checking for simple model fits + grid
# =====================================================

# Designed to be called in iPython notebook w/ inline plots
# Simply feed in your grid, supernova remnant, data, etc...

def check_eta2_simp(snr, kevs, data, eps, eta2_vals, mu_vals, fmt_vals):
    """Plot chi-squared space for eta2"""
    plt.figure(figsize=(6,4))
    for mu, fmt in zip(mu_vals, fmt_vals):
        b0_vals = []
        chisqr_vals = []
        for eta2 in eta2_vals:
            res = models.simple_fit(snr, kevs, data, eps, mu, eta2=eta2, eta2_free=False)
            b0_vals.append(res.params['B0'])
            chisqr_vals.append(res.chisqr)

        plt.plot(eta2_vals, chisqr_vals, fmt, lw=2)

    fplot(r'$\eta_2$ (-)', r'Best fit $\chi^2$', scales=('log','log'))
    plt.legend( tuple(r'$\mu = {:0.2f}$'.format(mu) for mu in mu_vals),
                loc='best')
    plt.show()

# ================================================
# Errors from confidence limits, simple model fits
# ================================================

# This is a simple reimplementation of `lmfit.conf_interval(...)`, using 
# a different metric for conf intervals. Compare results to be sure though

# get_ci_errors(...)
# get_ci_bounds(...) is a wrapper for get_bounds(...)
# one_dir_rootfind(...) is generalized, but convenient for chisqr thresholding


def get_ci_errors(ci_dict, par_name, ci_val=0.683):
    """Read out +/- errors from dict of CI bounds"""
    x_bnds = ci_dict[par_name]
    x_lo, x_hi = np.sort([tup[1] for tup in x_bnds if tup[0] == ci_val])
    x = [tup[1] for tup in x_bnds if tup[0] == 0.][0]
    return x - x_lo, x_hi - x


def get_ci_bounds(res, ci_vals=(0.683, 0.9)):
    """Manually compute confidence interval parameter values, for free
    parameters in res (lmfit.Minimizer)

    Input:
        res is the output from lmfit.minimize(...), ALREADY MINIMIZED
        f is the objective function being minimized
    Output:
        dict of ci bounds, structured like lmfit.conf_interval(...) output
    """
    ci2 = {}
    
    for pstr in [p for p in res.params if res.params[p].vary]:
        ci2[pstr] = [(0., res.params[pstr].value)]

        for ci_val in ci_vals:
            p_lo, p_hi = get_bounds(ci_val, res, pstr, eps=None, adapt=True,
                                    verbose=True)
            ci2[pstr].append((ci_val, p_lo))
            ci2[pstr].append((ci_val, p_hi))

        ci2[pstr].sort(key=lambda x: x[1])

    return ci2


def get_bounds(conf_intv, res, pstr, **kwargs):
    """Bounds for confidence interval conf_intv, within res.params[pstr] limits
    res must be already fitted, with fit parameters to vary in chisqr
    analysis already free.
    res will be modified in this process!
    
    Input: **kwargs are passed to one_dir_rootfind(...)
    Output: lower, upper bounds on parameter res.params[pstr]
    """

    # Rerun fit, store parameter values, stderrs, chisqr
    res.prepare_fit()
    res.leastsq()
    orig = lmfit.confidence.copy_vals(res.params)
    best_chisqr = res.chisqr

    # Fix parameter to be bounded throughout manual stepping
    res.params[pstr].vary = False  
    # Objective function for root finding.  Modifies res...
    # WARNING: other free parameters are constrained by their own limits!
    def chisqr_thresh(x):
        # Reset parameters to initial best fits, every time
        # seems like this reduces errors... but very ad hoc patch
        lmfit.confidence.restore_vals(orig, res.params)
        res.params[pstr].value = x
        res.prepare_fit()
        res.leastsq()
        return (res.chisqr - best_chisqr) - CI_2_CHISQR[conf_intv]

    # Find bounds if possible, using stored limits
    p_init = res.params[pstr].value
    p_lims = res.params[pstr].min, res.params[pstr].max
    p_lo = one_dir_rootfind(chisqr_thresh, p_init, p_lims[0], **kwargs)
    p_hi = one_dir_rootfind(chisqr_thresh, p_init, p_lims[1], **kwargs)

    # Restore values/chisqr in lmfit.Minimizer object
    lmfit.confidence.restore_vals(orig, res.params)
    res.params[pstr].vary = True
    res.chisqr = best_chisqr

    return p_lo, p_hi


def one_dir_rootfind(f, x_init, x_lim, eps=None, adapt=True, verbose=False):
    """Find root in function f from starting x with max/min limit
    Naive, brute force search in one direction, set by sgn(x_lim - x_init)
    If zero crossing is found, uses scipy.optimize.brentq to get root

    Designed specifically for chisqr bounds to compute confidence intervals

    Input: f has call signature f(x)
    Output: root within 2x machine precision, else x_lim if root was not found
    """
    if eps is None:
        eps = abs(0.01 * x_init) * np.sign(x_lim - x_init)
    else:
        eps = abs(eps) * np.sign(x_lim - x_init)

    prev_x = x_init
    prev_dist = f(x_init)
    start_sgn = np.sign(prev_dist)  # Allow crossings in either direction
    if start_sgn == 0:
        return x_init

    while np.sign(prev_dist) == start_sgn:
        # Break if we hit parameter limit without crossing threshold
        if np.abs(prev_x - x_lim) <= np.finfo(float).eps * 2:  # brentq rtol
            break

        x = prev_x + eps
        if (eps > 0 and x > x_lim) or (eps < 0 and x < x_lim):
            if verbose:
                print 'Hit parameter limit {}'.format(x_lim)
            x = x_lim

        dist = f(x)

        # Update step size if desired... arbitrary picks for adaptive step
        if adapt:
            if abs(dist - prev_dist) < 0.1:
                eps *= 1.5
            elif abs(dist - prev_dist) > 1:
                eps /= 1.5

        # Go back to start, where we can check if crossing was found or not
        prev_dist = dist
        prev_x = x

    # Broke out of loop (x = x_lim, prev_dist < 0),
    # or hit crossing right on (prev_dist == 0)
    if np.sign(prev_dist) == start_sgn:
        if verbose:
            print 'Crossing not found, hit limit at {}'.format(x)
        return x

    return sp.optimize.brentq(f, min(x_init,x), max(x_init,x))


if __name__ == '__main__':
    main()
