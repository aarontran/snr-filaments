"""
Code to execute model fits for measured FWHMs (using models, tables, etc
given in models_all.py

Here is where we should generate tables, figures, check effects of various
parameters. Plot chisquared space and check that grids are well sampled,
calculate errors from grids etc

Code should be SNR independent.
We want to be able to feed in any SNR, any data, and get the right fits.


Methods (for my own reference):

    TODO: fitting tables should spit to simple np.ndarray
    Array can be fed into 1. LaTeX converter, 2. iPython notebook display
    Currently output is just for iPython display.

    # Full model, fitting
    table_full(...)         generate tables/plots for full model, calls
    get_table_fit(...)      helpers to compute fits/errors
    get_table_err(...)

    # Full model, plotting
    check_eta2_grid(data, eps, tab, eta2, mus, fmts, inds=None)
    check_B0_grid(data, eps, tab, eta2, mus, fmts, inds=None)

    # Full model, grid (eta2dict) manipulation, get chisqr values
    grid_best(data, eps, eta2_dict, inds=None)
    grid_scan(data, eps, eta2_dict, inds=None)
    fwhm_scan(data, eps, fwhms, inds=None)

    # Simple model, fitting
    table_simp(...)         generate tables/plots for simple model
                            performs fits + error computation
    # Simple model, plotting
    check_eta2_simp(...)

    # Simple model, error bounds
    get_ci_errors(ci_dict, ...)
    get_ci_bounds(...)              wrapper for get_bounds
    get_bounds(...)
    one_dir_rootfind(...)           general, but meant for chisqr thresholding

Aaron Tran
2014 July 26
"""

from __future__ import division

#import cPickle as pickle
#from datetime import datetime
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import optimize
#import sys

from fplot import fplot
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


# TODO still in progress

def table_full(snr, kevs, data, eps, mus, fmts):
    """Make tables, main method.

    Find best grid value (could interp B0 here)
    Fit from best grid value, with multiple steps
    Anneal errors
    Spit out tables
    """
    table = []
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    for mu, fmt in zip(mus, fmts):

        # Find best grid value, interpolating in B0 if desired?
        eta2_vals, B0_vals, chisqr_vals = models.table_scan(data, eps, tab[mu], inds)
        ind = np.argmin(chisqr_vals)
        e2b = eta2_vals[ind], B0_vals[ind]

        # Fit from grid value
        # full_table_fit(...)

        # Get annealed errors (redundant call to table_scan)
        # TODO: SET RMINARC
        e2b_L, e2b_R = best_grid_pts(e2b[0], chisqr_vals[ind],
                                     snr, kevs, data, eps, mu, tab,
                                     conf_intv=0.683, inds=inds,
                                     fit_on_fly=True)

        # Assemble table for pretty-printing
        tr = ['{:0.2f}'.format(mu)]
        tr.append('{:0.3f} ({:0.3f}, {:0.3f})'.format(
                   e2b[0], e2b_L[0], e2b_R[0]))
        tr.append('{:0.2f} ({:0.2f}, {:0.2f})'.format(
                   e2b[1]*1e6, e2b_L[1]*1e6, e2b_R[1]*1e6))
        tr.append('{:0.4f}'.format(chisqr_vals[ind]))
        table.append(tr)

    return table, ax


def run_table_fit(snr, kevs, data, eps, mu_vals, fmt_vals, inds=None):
    """Execute table_fit for one set of data, w/ various mu values"""
    
    table = ListTable()
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    plt.figure(figsize=(6,4))
    plt.errorbar(kevs, data, yerr=eps, fmt='ok')

    for mu, fmt in zip(mu_vals, fmt_vals):

        print 'Starting fit from best grid point'
        # ALSO CHANGE RMINARC TO BE TIGHTER... BASED ON GRID FWHMS... may limit lmfit, that's okay.
        res, eta2_grid, B0_grid, chisqr_grid = models.table_fit(snr, kevs, data, eps, mu, tab,
                                                                inds=inds, epsfcn=1e-5)  # epsfcn!!!!!
        p = res.params
        print 'Best grid B0 = {:0.2f}, eta2 = {:0.2f} (chisqr = {:0.4f})'.format(B0_grid*1e6, eta2_grid, chisqr_grid)
        print 'Final fit chisqr = {:0.4f}'.format(res.chisqr)

        if chisqr_grid < res.chisqr:
            print 'Warning: lmfit worsened parameter fit.  Using best grid values, printed errors may be invalid'
            p['eta2'].value = eta2_grid
            p['B0'].value = B0_grid

        table_row = ['{:0.2f}'.format(mu)]
        table_row.append('{:0.2g} &plusmn; {:0.2g}'.format(p['eta2'].value, p['eta2'].stderr))
        table_row.append('{:0.3g} &plusmn; {:0.2g}'.format(p['B0'].value * 1e6, p['B0'].stderr * 1e6))
        table_row.append('{:0.4f}'.format(min(chisqr_grid, res.chisqr)))  # AD HOC...
        table.append(table_row)
        
        fwhms_m = models.width_cont(p, kevs, snr)  # Using actual model... can't overload w/ data points
        plt.plot(kevs, fwhms_m, fmt)

    display(table)
    fplot('Energy (keV)', 'FWHM (arcsec)', axargs='tight')
    plt.legend( tuple(['Data'] + [r'$\mu = {:0.2f}$'.format(mu) for mu in mu_vals]), loc='best')
    plt.show()



# ========================================================
# Grid fitting and parameter space checking for full model
# ========================================================

def get_table_fit(snr, kevs, data, eps, mu, tab, inds=None,
                   verbose=True, **lsq_kws):
    """Find best fit eta2, B0 for given mu value"""
    epsfcn_vals = [1e-10, 1e-8, 1e-6, 1e-4]  # MAGIC...

    # Initialize with grid values
    eta2, B0, fwhms, chisqr = models.grid_best(data, eps, tab[mu], inds=inds)
    eta2_best, B0_best = eta2, B0
    chisqr_best = chisqr

    for epsfcn in epsfcn_vals:
        if verbose:
            print 'Fitting from grid with epsfcn={:0.2g}'.format(epsfcn)

        res = models.full_fit(snr, kevs, data, eps, mu, eta2=eta2_grid, B0=B0_grid,
                              eta2_free=True, verbose=verbose,
                              epsfcn=epsfcn, **lsq_kws)

        if res.chisqr < chisqr_best:
            eta2_best = res.params['eta2'].value
            B0_best = res.params['B0'].value

    # Check if fitting improved grid value
    refit = True
    if chisqr_best == chisqr:
        if verbose:
            print 'lmfit did not improve best grid fit'
        refit = False

    return eta2_best, B0_best, chisqr_best, refit


def get_table_err(eta2_best, chisqr_min, snr, kevs, data, eps, mu, tab,
        inds=None, conf_intv=0.683, fit_on_fly=False):
        pass

# TODO still in progress

def best_grid_pts(eta2_best, chisqr_min, snr, kevs, data, eps, mu, tab, inds=None, conf_intv=0.683, fit_on_fly=False):
    """Estimate error in eta2/B0 simultaneously.  Threshold chisqr space in eta2 based on
    externally provided best fit eta2 + chisqr, and find most extreme points in grid bracketing
    this threshold.

    Do NOT naively trust function output.  Plot the data and ensure it looks
    sensible, as chi-squared space for full model fits is terribly behaved

    To get better estimates, at those points, perform a fit for the best value of B0.
    If that brings the point under the threshold, move on.  This will be rather slow.

    NOTE: don't be surprised if B0 from the "left" bound is larger than the best fit value...
    w/imperfect gridding, we don't get monotonicity in B0.
    """
    eta2_vals, B0_vals, chisqr_vals = models.table_scan(data, eps, tab[mu], inds)

    ind = np.searchsorted(eta2_vals, eta2_best)

    # Find where chisqr threshold for confidence interval is bracketed
    sgns = np.sign(chisqr_vals - (chisqr_min + CI_2_CHISQR[conf_intv]))
    crossings = np.where(np.diff(sgns))[0]  # First indices of crossings
    print 'minimum at {}'.format(ind)
    print 'crossings at {}'.format(crossings)

    flg = False  # DEBUGGING
    if any(crossings < ind):
        ind_left = crossings[0]
    else:
        ind_left = 0
        print 'bound not found (min eta2)'
        flg=True
    if any(crossings+1 > ind):
        ind_right = crossings[-1] + 1  # +1 to get right-side of crossing, + account for crossings AT best fit...
    else:
        ind_right = len(eta2_vals) - 1
        print 'bound not found (max eta2)'
        flg=True

    #if flg:
    #    print chisqr_vals - (chisqr_min + CI_2_CHISQR[conf_intv])
    #    print crossings

    e2b_left = eta2_vals[ind_left], B0_vals[ind_left]
    chisqr_left = chisqr_vals[ind_left]

    e2b_right = eta2_vals[ind_right], B0_vals[ind_right]
    chisqr_right = chisqr_vals[ind_right]

    print 'naive grid bounds (eta2): {:0.4f}, {:0.4f}'.format(e2b_left[0], e2b_right[0])
    print 'naive grid bounds (B0): {:0.4f}, {:0.4f}'.format(e2b_left[1]*1e6, e2b_right[1]*1e6)

    # Now, if there is leeway to explore...
    # ALSO CHANGE RMINARC TO BE TIGHTER...
    if fit_on_fly:
        if ind_left > 0:
            print 'Checking left error bound (original bound chisqr = {}, B0 = {:0.4f})'.format(chisqr_left, e2b_left[1]*1e6)

            # Need to check results are independent of epsfcn?!
            # Should have step size comparable to, or just smaller than, grid spacing in B0 -- to ensure
            # we are just exploring that in-between-points area of parameter space.  This we can do.
            res = models.full_fit(snr, kevs, data, eps, mu, eta2=e2b_left[0], B0=e2b_left[1], epsfcn=1e-6)
            while ind_left > 0 and res.chisqr - (chisqr_min + CI_2_CHISQR[conf_intv]) < 0:
                print 'Threshold moved by one (new best B0 here = {:0.4f}) chisqr was {}; min is {}'.format(res.params['B0'].value*1e6,
                                                                                                       res.chisqr, chisqr_min)
                ind_left -= 1
                e2b_left = eta2_vals[ind_left], B0_vals[ind_left]
                res = models.full_fit(snr, kevs, data, eps, mu, eta2=e2b_left[0], B0=e2b_left[1], epsfcn=1e-6)
            e2b_left = (e2b_left[0], res.params['B0'].value)

        if ind_right < len(eta2_vals) - 1:
            print 'Checking right error bound (original bound chisqr = {}, B0 = {:0.4f})'.format(chisqr_right, e2b_right[1]*1e6)
            res = models.full_fit(snr, kevs, data, eps, mu, eta2=e2b_right[0], B0=e2b_right[1], epsfcn=1e-6)
            while ind_right < len(eta2_vals) - 1 and res.chisqr - (chisqr_min + CI_2_CHISQR[conf_intv]) < 0:
                print 'Threshold moved by one (new best B0 here = {:0.4f}) chisqr was {}; min is {}'.format(res.params['B0'].value*1e6,
                                                                                                       res.chisqr, chisqr_min)
                ind_right += 1
                e2b_right = eta2_vals[ind_right], B0_vals[ind_right]
                res = models.full_fit(snr, kevs, data, eps, mu, eta2=e2b_right[0], B0=e2b_right[1], epsfcn=1e-6)
            e2b_right = (e2b_right[0], res.params['B0'].value)

    print 'final grid bounds (eta2): {:0.4f}, {:0.4f}'.format(e2b_left[0], e2b_right[0])
    print 'final grid bounds (B0): {:0.4f}, {:0.4f}'.format(e2b_left[1]*1e6, e2b_right[1]*1e6)

    return e2b_left, e2b_right



# ======================================================
# Parameter space checking for full model grid (no fits)
# ======================================================

def check_eta2_grid(data, eps, tab, mus, fmts, inds=None):
    """Plot chi-squared space for eta2 based on grid values
    JACKSON POLLOCK DOES ASTROPHYSICS
    """
    ax = plt.gca()
    for mu, fmt in zip(mus, fmts):
        eta2s, B0s, fwhms, chisqrs = grid_scan(data, eps, tab[mu], inds)
        ax.plot(eta2s, chisqrs, fmt)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(tuple(r'$\mu = {:0.2f}$'.format(mu) for mu in mus), loc='best')
    return ax


def check_B0_grid(data, eps, tab, eta2, mus, fmts, inds=None):
    """Plot chisqr as a function of B0 for arbitrary (mu, eta2), with
    multiple mu values on a single plot.
    Verify that we can resolve chisq minimum in B0
    """
    ax = plt.gca()
    for mu, fmt in zip(mus, fmts):
        B0_vals, fwhms = tab[mu][eta2]
        chisqr_vals = fwhm_scan(data, eps, fwhms, inds=inds)
        ax.plot(B0_vals*1e6, chisqr_vals, fmt)
    return ax


# ===========================================================
# Grid manipulation -- compute chisqrs and find best grid pts
# ===========================================================

def grid_best(data, eps, eta2_dict, inds=None):
    """Find the best fit to data/eps in grid specified by eta2_dict
    A silly wrapper... but, sometimes we want grid_scan output
    directly for plots.
    """
    eta2s, B0s, fwhms, chisqrs = grid_scan(data, eps, eta2_dict, inds=inds)
    ind = np.argmin(chisqrs)
    return eta2s[ind], B0s[ind], fwhms[ind], chisqrs[ind]


def grid_scan(data, eps, eta2_dict, inds=None):
    """Compute chisqr over grid for provided data/eps, subsetting with inds
    if necessary (bands in data, eps should match grid)

    Input: data, eps must be np.ndarray
    Output: SORTED eta2 values with corresponding best B0, fwhms, chisqr
    """
    eta2_vals = np.sort(eta2_dict.keys())
    B0_vals = np.array([])
    fwhm_vals = np.array([])
    chisqr_vals = np.array([])

    for eta2 in eta2_vals:
        B0_pts, fwhms = eta2_dict[eta2]

        chisqr_pts = fwhm_scan(data, eps, fwhms, inds=inds)

        ind_min = np.argmin(chisqr_pts)  # Could interpolate a better B0 fit
        B0_vals = np.append(B0_vals, B0_pts[ind_min])
        fwhm_vals = np.append(fwhm_vals, fwhms[ind_min])
        chisqr_vals = np.append(chisqr_vals, chisqr_pts[ind_min])

    return eta2_vals, B0_vals, fwhm_vals, chisqr_vals


def fwhm_scan(data, eps, fwhms, inds=None):
    """Compute chisqr values over B0 mesh"""
    if inds:  # Subset if need
        data, eps = data[inds], eps[inds]
        fwhms = fwhms[:,inds]

    return map(lambda x: models.chi_squared(data, eps, x), fwhms)



# =====================================================
# Parameter space checking for simple model fits + grid
# =====================================================

# Designed to be called in iPython notebook w/ inline plots
# Simply feed in your grid, supernova remnant, data, etc...

def table_simp(snr, kevs, data, eps, mus, fmts):
    """Generate Table 7 of Sean's paper with better 1-sigma errors + plot"""

    table = []  # To be used with ListTable in iPython notebook output
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    ax = plt.gca()
    ax.errorbar(kevs, data, yerr=eps, fmt='ok')

    for mu, fmt in zip(mus, fmts):
        # Actual computations / fitting
        res = models.simple_fit(snr, kevs, data, eps, mu)
        p = res.params
        ci2 = get_ci_bounds(res, ci_vals=(0.683,))
        p['eta2'].better_err = get_ci_errors(ci2, 'eta2', ci_val=0.683)
        p['B0'].better_err = get_ci_errors(ci2, 'B0', ci_val=0.683)

        tr = ['{:0.2f}'.format(mu)]
        tr.append('{:0.2g} &plusmn; {:0.2g} (or, +{:0.2g}/-{:0.2g})'.format(
                  p['eta2'].value, p['eta2'].stderr,
                  p['eta2'].better_err[1], p['eta2'].better_err[0]))
        tr.append('{:0.3g} &plusmn; {:0.2g} (or, +{:0.2g}/-{:0.2g})'.format(
                  p['B0'].value * 1e6, p['B0'].stderr * 1e6,
                  p['B0'].better_err[1]*1e6, p['B0'].better_err[0]*1e6))
        tr.append('{:0.4f}'.format(res.chisqr))
        table.append(tr)

        kevs_m = np.linspace(kevs[0]-0.2, kevs[-1]+0.2, 100)
        fwhms_m = models.width_dump(p, kevs_m, snr)
        ax.plot(kevs_m, fwhms_m, fmt)

    fplot('Energy (keV)', 'FWHM (arcsec)', axargs='tight')
    ax.legend(tuple(['Data'] + [r'$\mu = {:0.2f}$'.format(mu)
                                 for mu in mus]), loc='best')
    return table, ax


def check_eta2_simp(snr, kevs, data, eps, eta2_vals, mus, fmts, inds=None):
    """Plot chi-squared space for eta2"""
    ax = plt.gca()

    if inds is not None:
        kevs, data, eps = kevs[inds], data[inds], eps[inds]

    for mu, fmt in zip(mus, fmts):
        b0_vals = []
        chisqr_vals = []
        for eta2 in eta2_vals:
            res = models.simple_fit(snr, kevs, data, eps, mu,
                                    eta2=eta2, eta2_free=False)
            b0_vals.append(res.params['B0'])
            chisqr_vals.append(res.chisqr)

        ax.plot(eta2_vals, chisqr_vals, fmt, lw=2)

    fplot(r'$\eta_2$ (-)', r'Best fit $\chi^2$', scales=('log','log'))
    ax.legend( tuple(r'$\mu = {:0.2f}$'.format(mu) for mu in mus),
                loc='best')
    return ax


# ================================================
# Errors from confidence limits, simple model fits
# ================================================

# Simple reimplementation of `lmfit.conf_interval(...)`, using a
# different metric for conf intervals. Compare results to be sure though


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
