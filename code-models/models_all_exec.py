"""
Code to execute model fits for measured FWHMs

Compute best-fit parameters and errors
Generate tables, figures, check effects of various parameters.

Aaron Tran
2014 July 26
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import optimize as spopt

import lmfit

from fplot import fplot
from np2latex import LatexTable
import models_all as models

# TODO not sure of best way to store/hold this
# Table from Press et al., Section 14.5 (1st ed.),
# 15.6 (3rd ed.); see also Avni, Y. (1976, ApJ)
# Would be good to compute properly in a function
CI_2_CHISQR = {}
CI_2_CHISQR[0.683] = 1.0
CI_2_CHISQR[0.90] = 2.71
CI_2_CHISQR[0.9545] = 4.0
CI_2_CHISQR[0.99] = 6.63


def main():
    """For debugging/testing stuff.  Not currently in use"""
    pass

#TODO: fitting methods should spit to simple np.ndarray, instead of wrangling
#      with all of this display / plotting stuff
#
#      Input: snr, kevs, data+eps, tab, mus
#      Output: np.ndarray of [mu, eta2, eta2 errors, B0, B0 errors, chisqr]
#
#      Then, write methods (another module) to generate
#      LaTeX, iPython notebook tables from np.ndarrays

# ===========================
# Full model fitting/plotting
# ===========================

def table_full(snr, kevs, data, eps, tab, mus, fmts, ax=None, inds=None,
               verbose=True):
    """Perform full model fits, make tables and plots of outputs

    Find best grid value (could interp B0 here)
    Fit from best grid value, with multiple steps
    Anneal errors
    Spit out tables
    """
    table = []
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    ltab = LatexTable([r'$\mu$ (-)', r'$\eta_2$ (-)', r'$B_0$ ($\mu$G)',
                       r'$\chi^2$'],
                      ['{:0.2f}', 2, 2, '{:0.4f}'],
                       'LaTeX table')  # BAD...

    if ax is None:
        ax = plt.gca()
    ax.errorbar(kevs, data, yerr=eps, fmt='ok')

    for mu, fmt in zip(mus, fmts):
        # Fit from best grid value NOTE FITTING SWITCH
        eta2, B0, chisqr = get_table_fit(snr, kevs, data, eps, mu, tab,
                                         inds=inds, verbose=verbose)
        #print 'only using best grid value'
        #eta2, B0, asdfjkl_dont_care, chisqr = grid_best(data,eps,tab[mu],inds)

        if verbose:
            print 'Best fit values: eta2={:0.4f}, B0={:0.4f}'.format(eta2,B0*1e6)

        # Get annealed errors
        print 'computing errors'
        ci = get_table_err(eta2, B0, chisqr, snr, kevs, data, eps, mu, tab,
                           conf_intv=0.683, inds=inds, verbose=verbose,
                           anneal=True)  # NOTE ANNEALING SWITCH
        err_eta2 = get_ci_errors(ci, 'eta2', ci_val=0.683)
        err_B0 = get_ci_errors(ci, 'B0', ci_val=0.683)

        # Assemble table for pretty-printing
        tr = ['{:0.2f}'.format(mu)]
        tr.append('{:0.3f} (+{:0.3f}/-{:0.3f})'.format(
                   eta2, err_eta2[1], err_eta2[0]))
        tr.append('{:0.2f} (+{:0.2f}/-{:0.2f})'.format(
                   B0*1e6, err_B0[1]*1e6, err_B0[0]*1e6))
        tr.append('{:0.4f}'.format(chisqr))
        table.append(tr)

        # Assemble LaTeX table
        ltab.add_row(mu, eta2, err_eta2[1], err_eta2[0],
                     B0*1e6, err_B0[1]*1e6, err_B0[0]*1e6,
                     chisqr)

        # Plot
        p = lmfit.Parameters()
        p.add('B0',value=B0)
        p.add('eta2',value=eta2)
        p.add('mu',value=mu)
        fwhms_m = models.width_cont(p, kevs, snr)  # Here's a function call
        print fwhms_m
        ax.plot(kevs, fwhms_m, fmt)
        ax.legend(tuple(['Data'] + [r'$\mu = {:0.2f}$'.format(mu) for
                                    mu in mus]), loc='best')
        fplot('Energy (keV)', 'FWHM (arcsec)', ax=ax, axargs='tight')

    return table, ax, ltab

# ========================================================
# Grid fitting and parameter space checking for full model
# ========================================================

def get_table_fit(snr, kevs, data, eps, mu, tab, inds=None,
                   verbose=True, **lsq_kws):
    """Find best fit eta2, B0 for given mu value (need mu for fit)"""
    f_rmin = 1.2  # TODO MAGIC...

    # Initialize with best grid values
    eta2, B0, fwhms, chisqr = grid_best(data, eps, tab[mu], inds=inds)
    if verbose:
        print 'chisqr from tabulated fwhms: {}'.format(chisqr)

    # methods for scanning grid incorporate inds.  But, full_fit/width_cont
    # do not do so. (THIS MUST COME AFTER grid_best(...) call...
    # TODO BAD DESIGN)
    if inds is not None:
        kevs = kevs[inds]
        data, eps = data[inds], eps[inds]

    # Update fwhms, chisqr as table values may be inaccurate/out-of-date
    fwhms = models.full_width(snr, kevs, mu, eta2, B0, verbose=False)
    chisqr = models.chi_squared(data, eps, fwhms)
    if verbose:
        print 'chisqr w/ more accurate fwhms: {}'.format(chisqr)

    # Now, fit from best grid point
    # TODO: save lmfit standard error estimates, for comparison...
    res = models.full_fit(snr, kevs, data, eps, mu, eta2=eta2, B0=B0,
                          eta2_free=True, verbose=verbose,
                          #rminarc=fwhms*f_rmin,
                          epsfcn=1e-6, **lsq_kws)  # TODO epsfcn magic number
    if verbose:
        print 'new chisqr = {}, grid chisqr = {}'.format(res.chisqr, chisqr)
    if res.chisqr < chisqr:
        eta2 = res.params['eta2'].value
        B0 = res.params['B0'].value
        chisqr = res.chisqr
        print 'Found better fit at eta2={:0.4f}, B0={:0.4f}'.format(eta2,B0*1e6)
    else:
        if verbose:
            print 'Error: failed to improve best grid values (?!)'

    return eta2, B0, chisqr

# TODO ASSUMPTION -- extreme values of eta2 also give correct B0 bounds
# But this has a good chance of being incorrect.  Should repeat this process
# for B0... argh :(

# NOTE currently, spopt in between "bracketing" values FAILS because
# pre-cached FWHMs are old and give inaccurate FWHMs.  Thus, grid bracketing
# values may not bracket the threshold when tables are recalculated.
# (this isn't an issue if you anneal more than one point)
# I won't fix this yet -- will redesign code first.

def get_table_err(eta2, B0, chisqr_min, snr, kevs, data, eps, mu, tab,
                  inds=None, verbose=True, anneal=True, conf_intv=0.683):
    """Estimate error in eta2, B0 from pre-tabulated grid.

    Threshold eta2-chisqr space based on given best fit eta2 + chisqr, and
    find most extreme points in grid bracketing this threshold.
    "Anneal" errors by fitting for best B0.  If this drops chisqr below
    threshold, move to next eta2 (on left or right)

    mu, tab passed in separately (instead of as eta2_dict) because
    we need access to mu value for full model calculations/fitting

    Input:
    Output:
    """
    # Ad hoc patch (stackoverflow.com/a/5980173)
    if verbose:
        def vprint(*args):
            for arg in args:
                print arg,
            print
    else:
        vprint = lambda *a: None

    # Load grid, locate search start point
    eta2grid, B0grid, fwhmsgrid, chisqrgrid = grid_scan(data,eps,tab[mu],inds)
    idx = np.searchsorted(eta2grid, eta2)  # Location of best eta2
    def chisqr_thresh(x):
        return x - (chisqr_min + CI_2_CHISQR[conf_intv])

    # Find where conf_intv chisqr threshold is bracketed in eta2 grid
    # If crossings are not found, start from grid edges
    sgns = np.sign(chisqr_thresh(chisqrgrid))
    crossings = np.where(np.diff(sgns))[0]  # Left indices of crossings
    idx_left = crossings[0] if any(crossings < idx) else 0
    # crossings[-1]+1 gets right index of crossing
    idx_right = crossings[-1]+1 if any(crossings+1 > idx) else len(eta2grid)-1

    # Check (anneal) grid chisqrs and find bracketing (eta2,B0) pts, then
    # use one_dir_root(...) or spopt.brentq(...) to get more precise limit
    eta2bnd = (0., 1e5)  # TODO ARBITRARY/MAGIC LIMITS on eta2 search

    class Fitter(object):  # Define here to access snr/kevs/etc
        """Stores prev B0 value for next thresh fit; inits with best fit B0"""
        def __init__(self):
            self.B0_prev = B0
        def eta2_thresh(self, eta2_test):
            # TODO set rminarc, epsfcn, iradmax, resolution etc appropriately
            # rminarc = fwhmsgrid[idx] * f_rmin
            res = models.full_fit(snr, kevs, data, eps, mu, eta2=eta2_test,
                                  B0=self.B0_prev, eta2_free=False,
                                  verbose=verbose)
            self.B0_prev = res.params['B0'].value  # Store previous value
            return chisqr_thresh(res.chisqr)

    if anneal:  # TODO err_from_grid design is bad, just checking that it works
        eta2_L, B0_L = err_from_grid(idx_left, -1, Fitter(), vprint,
                                     eta2grid, eta2bnd[0])
        eta2_R, B0_R = err_from_grid(idx_right, +1, Fitter(), vprint,
                                     eta2grid, eta2bnd[1])
    else:
        vprint('Using naive grid bounds (no annealing)')
        eta2_L, B0_L = eta2grid[idx_left], B0grid[idx_left]
        eta2_R, B0_R = eta2grid[idx_right], B0grid[idx_right]
        
    # Output in format analogous to lmfit.conf_interval
    ci = {}
    ci['eta2'] = [(conf_intv, eta2_L), (conf_intv, eta2_R), (0., eta2)]
    ci['B0'] = [(conf_intv, B0_L), (conf_intv, B0_R), (0., B0)]
    return ci

def err_from_grid(idx, sgn, fobj, vprint, eta2grid, eta2_search_bound):
    """Get error"""
    sgn = np.sign(sgn)
    if sgn < 0:
        vprint('Annealing on left')
    else:
        vprint('Annealing on right')
    idx = anneal_grid(idx, sgn, eta2grid, fobj.eta2_thresh)
    vprint('Done annealing')
    if idx is None:  # Hit grid edge
        edge = 0 if sgn<0 else -1
        vprint('Hit grid edge {}'.format(eta2grid[edge]))
        eta2_lim = one_dir_root(fobj.eta2_thresh, eta2grid[edge],
                                eta2_search_bound, verbose=True)
        if eta2_lim == eta2_search_bound:
            print 'Warning: hit search bound {}'.format(eta2_lim)
    else:
        vprint('Found bracketing values')
        # On left (sgn=-1) we want idx, idx+1
        # on right(sgn=+1) we want idx-1, idx
        eta2_lim = spopt.brentq(fobj.eta2_thresh,
                                eta2grid[min(idx, idx - sgn)],
                                eta2grid[max(idx, idx - sgn)])
    return eta2_lim, fobj.B0_prev

def anneal_grid(ind, sgn, eta2_vals, f_thresh):
    """Find grid index just above threshold; return None if not found"""
    f_res = f_thresh(eta2_vals[ind])  # Check first index, ind
    while f_res <= 0 and (0 < ind < len(eta2_vals)-1):
        print 'Moved threshold by one'
        ind += np.sign(sgn)
        f_res = f_thresh(eta2_vals[ind])
    return ind if f_res > 0 else None

# ===========================================================
# Grid manipulation -- compute chisqrs and find best grid pts
# ===========================================================

def grid_best(data, eps, eta2_dict, inds=None):
    """Find best fit to data/eps in (eta2, B0) grid of pretabulated FWHMs

    Input:
        data, eps (np.ndarray): data and errors to fit
        inds (np.ndarray): use only a subset of energy bands if desired
    Output:
        eta2 (float), B0 (float), fwhms (np.ndarray), chisqr (float)
    """
    eta2s, B0s, fwhms, chisqrs = grid_scan(data, eps, eta2_dict, inds)
    ind = np.argmin(chisqrs)
    return eta2s[ind], B0s[ind], fwhms[ind], chisqrs[ind]

def grid_scan(data, eps, eta2_dict, inds=None):
    """Find best-fit B0s, FWHMs, chisqrs for eta2 grid

    Compute chisqr over (eta2, B0) grid for provided data/eps,
    subset with inds if necessary.  Cull best B0 value for each eta2.

    Input: data, eps (np.ndarray); eta2_dict (dict) grid of eta2, B0, FWHMs
    Output: SORTED eta2 values with best B0, fwhms, chisqr
    """
    eta2_vals = np.sort(eta2_dict.keys())
    B0_vals = np.array([])
    fwhm_vals = []  # Easier to build up, but not efficient
    chisqr_vals = np.array([])

    for eta2 in eta2_vals:
        B0_pts, fwhms = eta2_dict[eta2]
        chisqr_pts = fwhm_scan(data, eps, fwhms, inds)

        # Previously tried interpolating better B0 fit, now scrapped
        ind_min = np.argmin(chisqr_pts)
        B0_vals = np.append(B0_vals, B0_pts[ind_min])
        fwhm_vals.append(fwhms[ind_min])
        chisqr_vals = np.append(chisqr_vals, chisqr_pts[ind_min])

    return eta2_vals, B0_vals, np.array(fwhm_vals), chisqr_vals

def fwhm_scan(data, eps, fwhms, inds=None):
    """Compute chisqr values over pre-computed FWHMs (from B0 mesh)"""
    if inds is not None:  # Subset if needed
        data, eps = data[inds], eps[inds]
        fwhms = fwhms[:,inds]

    return map(lambda x: models.chi_squared(data, eps, x), fwhms)


# ========================
# Fitting for simple model
# ========================

def table_simp(snr, kevs, data, eps, mus, fmts, inds=None):
    """Generate Table 7 of Sean's paper with better 1-sigma errors + plot"""

    table = []  # To be used with ListTable in iPython notebook output
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    ltab = LatexTable([r'$\mu$ (-)', r'$\eta_2$ (-)', r'$B_0$ ($\mu$G)',
                       r'$\chi^2$'],
                      ['{:0.2f}', 1, 1, '{:0.4f}'],
                       'LaTeX table')  # BAD...

    ax = plt.gca()
    ax.errorbar(kevs, data, yerr=eps, fmt='ok')

    if inds is not None:
        kevs, data, eps = kevs[inds], data[inds], eps[inds]

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

        ltab.add_row(mu, p['eta2'].value, p['eta2'].stderr,
                     p['B0'].value * 1e6, p['B0'].stderr * 1e6,
                     res.chisqr)

        kevs_m = np.linspace(kevs[0]-0.2, kevs[-1]+0.2, 100)
        fwhms_m = models.width_dump(p, kevs_m, snr)
        ax.plot(kevs_m, fwhms_m, fmt)

    fplot('Energy (keV)', 'FWHM (arcsec)', axargs='tight')
    ax.legend(tuple(['Data'] + [r'$\mu = {:0.2f}$'.format(mu)
                                 for mu in mus]), loc='best')
    return table, ax, ltab

# =========================
# Errors, simple model fits
# =========================
# Similar to `lmfit.conf_interval(...)`, but different metric for conf intvs

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

    Input: **kwargs are passed to one_dir_root(...)
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
    p_lo = one_dir_root(chisqr_thresh, p_init, p_lims[0], **kwargs)
    p_hi = one_dir_root(chisqr_thresh, p_init, p_lims[1], **kwargs)

    # Restore values/chisqr in lmfit.Minimizer object
    lmfit.confidence.restore_vals(orig, res.params)
    res.params[pstr].vary = True
    res.chisqr = best_chisqr

    return p_lo, p_hi

# TODO can be optimized, look inside this method if folded
def one_dir_root(f, x_init, x_lim, eps=None, adapt=True, verbose=False):
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

    # TODO instead of searching between x and x_init, search btwn x and prev x,
    # which should also bracket threshold.  Slightly more efficient/faster
    # (especially if using one_dir_root for full model error analysis)
    return sp.optimize.brentq(f, min(x_init,x), max(x_init,x))

# ========================
# Parameter space checking
# ========================

def check_eta2_grid(data, eps, tab, mus, fmts, inds=None):
    """Plot chisqr-eta2 space from (eta2,B0) grids for multiple mu values
    JACKSON POLLOCK DOES ASTROPHYSICS

    Input:
        data, eps, tab: as usual
        mus: list of mu values to check/plot
        fmts: list of linespec strings for plots
        inds: subset if desired
    Output:
        matplotlib.Axis object
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
    """Plot chisqr(B0) from grid for fixed eta2, multiple mu values

    Verify that we can resolve chisq minimum in B0 grid
    
    Input:
        data, eps, tab: as usual
        eta2: where in grid shall we slice?
        mus: list of mu values to check/plot
        fmts: list of linespec strings for plots
        inds: subset if desired
    Output:
        matplotlib.Axis object
    """
    ax = plt.gca()
    for mu, fmt in zip(mus, fmts):
        B0_vals, fwhms = tab[mu][eta2]
        chisqr_vals = fwhm_scan(data, eps, fwhms, inds=inds)
        ax.plot(B0_vals*1e6, chisqr_vals, fmt)
    return ax

def check_eta2_simp(snr, kevs, data, eps, eta2_vals, mus, fmts, inds=None):
    """Plot chi-squared space for eta2, based on simple model"""
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


if __name__ == '__main__':
    main()
