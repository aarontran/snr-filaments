"""
Code to execute model fits for measured FWHMs
Compute best-fit parameters and errors

Responsibility for various kwargs split between here and models.py

Aaron Tran
August 2014
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import optimize as spopt

import lmfit

import models

from fplot import fplot

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
    pass

class Fitter(object):

    def __init__(self, snr, kevs, data, eps, tab, inds=None, verbose=True):
        """Fitter w/SNR and measured data to run full/simple model fits.

        Use setters to update any of kevs/data/eps (does indexing for you)
        Updating kevs requires new indices as well
        """
        self.snr = snr
        self.tab = tab

        self._okevs = kevs
        self._odata = data
        self._oeps = eps

        self.set_inds(inds)
        self.set_verbose(verbose)
    # self.snr, self.inds, self.kevs, self.data, self.eps
    # self._verbose, self._vprint
    # Fitter doesn't keep a record of previous fits

    # -------
    # Setters
    # -------

    # Could improve these using hooks (similar to lmfit parameters.py)
    # but, really not important

    def set_inds(self, inds=None):
        """Update energy band indices to be used.  None resets to all"""
        self.inds = inds
        if inds is not None:
            self.kevs = self._okevs[inds]
            self.data = self._odata[inds]
            self.eps = self._oeps[inds]
        else:
            self.kevs = self._okevs
            self.data = self._odata
            self.eps = self._oeps

    def set_data_err(self, data, eps):
        """Update data/eps to be fitted, respecting inds"""
        self._odata = data
        self._oeps = eps
        self.set_inds(self.inds)

    def set_kevs(self, kevs, inds=None):
        """Update kevs; resets inds by default"""
        self._okevs = kevs
        self.set_inds(inds)

    def set_verbose(self, verbose):
        """Set verbosity (_verbose, _vprint)"""
        self._verbose = verbose
        if verbose:
            def _vprint(*args):
                for arg in args:
                    print arg,
                print
            self._vprint = _vprint
        else:
            self._vprint = lambda *a: None

    # -------------------------------------
    # Wrappers for quick width computations
    # -------------------------------------
    # Useful for one-off debugging, interactive work, etc.

    def chisqr_full(self, mu, eta2, B0, **kwargs):
        """chisqr for specified full model, for stored SNR and FWHM data

        **kwargs passed directly to full model code. kws (as of 2014 Sep 12):
            verbose, rminarc, icut, irmax, iradmax, ixmax,
            irad_adapt, irad_adapt_f
        """
        fwhm_m = width_full(mu, eta2, B0, **kwargs)

        return models.chi_squared(self.data, self.eps, fwhm_m)

    def chisqr_simp(self, mu, eta2, B0):
        """chisqr for specified simple model, for stored SNR/FWHM data"""
        fwhm_m = width_simp(mu, eta2, B0)
        return models.chi_squared(self.data, self.eps, fwhm_m)

    def width_full(self, mu, eta2, B0, **kwargs):
        """Full model widths for stored SNR and energy bands
        Wrapper for models.width_cont, can pass through intensities/etc as well

        **kwargs passed directly to full model code (see models.width_cont)
        """
        return models.full_width(self.snr, self.kevs, mu, eta2, B0, **kwargs)

    def width_simp(self, mu, eta2, B0):
        """Simple model widths for stored SNR and energy bands"""
        return models.simple_width(self.snr, self.kevs, mu, eta2, B0)

    # -----------------------------
    # Wrappers for fitting routines
    # -----------------------------

    def get_fit(self, mu, appr, **kwargs):
        """Simple or full fit + manual 1-sigma errors
        **kwargs passed to models.full_fit, if appr='full'
        In particular, kws include: model_kws, method, epsfcn, maxfev, etc.

        Here, I specifically disable redchi-scaling for the covar matrix
        unless specified otherwise (use kwarg 'scale_covar')
        """
        if 'scale_covar' not in kwargs:
            kwargs['scale_covar'] = False

        self._vprint('Best {} fit, mu = {}'.format(appr, mu))

        if appr == 'simp':
            res = models.simple_fit(self.snr, self.kevs, self.data, self.eps,
                                    mu, **kwargs)

        elif appr == 'full':

            if ('eta2' not in kwargs) or ('B0' not in kwargs):

                eta2_vals, B0_vals, fwhm_vals, chisqr_vals = self.grid_scan(mu)

                if 'eta2' in kwargs:
                    self._vprint('Starting from provided eta2, best grid B0')
                    ind = np.searchsorted(eta2_vals, kwargs['eta2'])
                    if ind >= len(eta2_vals):  # If eta2 outside gridded values
                        ind = len(eta2_vals) - 1  # guess B0 fr largest eta2
                    kwargs['B0'] = B0_vals[ind]
                elif 'B0' in kwargs:
                    self._vprint('Starting from provided B0, best grid eta2')
                    msk = np.argsort(B0_vals)
                    eta2_vals = eta2_vals[msk]
                    B0_vals = B0_vals[msk]
                    ind = np.searchsorted(B0_vals, kwargs['B0'])
                    if ind >= len(B0_vals):  # As above, if B0 outside all best
                        ind = len(B0_vals) - 1  # values, get the closest one
                    kwargs['eta2'] = eta2_vals[ind]
                else:
                    self._vprint('Starting from best grid eta2, B0')
                    ind = np.argmin(chisqr_vals)
                    kwargs['eta2'] = eta2_vals[ind]
                    kwargs['B0'] = B0_vals[ind]

            else:
                self._vprint('Starting from user provided eta2, B0')

            self._vprint('Initial guesses:',
                         'eta2 = {:0.3f},'.format(kwargs['eta2']),
                         'B0 = {:0.3f}'.format(kwargs['B0']*1e6))
            res = models.full_fit(self.snr, self.kevs, self.data, self.eps,
                                  mu, **kwargs)

        else:
            raise Exception('Invalid model choice')

        self._vprint('Done fitting')

        return res

    def get_errs(self, res, appr, method='manual', ci_vals=(0.683,), **kwargs):
        """Find dictionary of confidence intervals

        Input:
            res = lmfit.Minimizer() from fit
            appr = 'full' or 'simp'
            method = 'manual' (default), 'lmfit', or 'stderr'

            **kwargs are passed to relevant functions
                method='lmfit': maxiter=200, prob_func=None
                                (see lmfit.conf_interval docstring)
                method='manual': adapt=True, anneal=True, eps=None;
                                 remaining **kwargs sent to brentq
                                 (anneal keyword only if appr='full')
        Output:
            confidence interval dictionary
        """
        self._vprint('Finding {} fit errors; method={}'.format(appr, method))

        if method == 'lmfit':
            return lmfit.conf_interval(res, sigmas=ci_vals,
                                       verbose=self._verbose, **kwargs)
        if method == 'stderr':
            def f(conf_intv, res, pstr):
                if conf_intv != 0.683:
                    print ('Warning: requested CI={}. But, stderr only '
                           'gives 68.3% errors').format(conf_intv)
                v, s = res.params[pstr].value, res.params[pstr].stderr
                return v - s, v + s
        elif method == 'manual':
            f = lambda *args: self.get_bounds(*args, appr=appr,
                                verbose=self._verbose, **kwargs)
        else:
            raise Exception('ERROR: method must be one of lmfit/stderr/manual')

        ci = build_ci_dict(res, f, ci_vals=ci_vals)
        if self._verbose:
            self._vprint('Finished finding fit errors:')
            print lmfit.printfuncs.ci_report(ci)
            self._vprint('')
        return ci

    # ----------------------------
    # Model specific error finding
    # ----------------------------

    def get_bounds(self, conf_intv, res, pstr, appr='full', anneal=True,
                   **kwargs):
        """Estimate error in eta2, B0 from pre-tabulated grid.

        Threshold eta2-chisqr space based on given best fit eta2 + chisqr, and
        find most extreme points in grid bracketing this threshold.
        "Anneal" errors by fitting for best B0.  If this drops chisqr below
        threshold, move to next eta2 (on left or right)

        Input: res is best fit lmfit.Minimizer
            **kwargs passed through to one_dir_root (and brentq)
        Output: limits on parameter pstr
        """
        orig = lmfit.confidence.copy_vals(res.params)
        best_chisqr = res.chisqr

        res.params[pstr].vary = False
        def chisqr_thresh(x):
            return x - (best_chisqr + CI_2_CHISQR[conf_intv])
        # Objective function for root-finding.  Modifies res, doesn't restore
        # so, when stepping on grid / root-finding, it keeps memory
        # of its previous fit state
        def f_thresh(x):  # TODO other free params constrained by limits!
            res.params[pstr].value = x
            res.prepare_fit(pstr)
            res.leastsq()
            return chisqr_thresh(res.chisqr)

        # Best fit quantities
        p_init = res.params[pstr].value
        p_min, p_max = res.params[pstr].min, res.params[pstr].max

        self._vprint('Finding error bounds on {} = {}'.format(pstr, p_init))

        if appr == 'simp':
            self._vprint('Bounding below, limit = {}'.format(p_min))
            p_lo = one_dir_root(f_thresh, p_init, p_min, **kwargs)

            self._vprint('Bounding above, limit = {}'.format(p_max))
            lmfit.confidence.restore_vals(orig, res.params)
            p_hi = one_dir_root(f_thresh, p_init, p_max, **kwargs)

        elif appr == 'full':
            # Get xgrid from self.tab + lo/hi indices to search from
            _stuff = self.prep_grid_anneal(chisqr_thresh, res, pstr)
            xgrid, idx_lo, idx_hi, others_lo, others_hi = _stuff  # yeah...

            if anneal:
                self._vprint('Bounding below, limit = {}'.format(p_min))
                for opstr, op in others_lo.items():  # Load other params
                    res.params[opstr].value = op
                p_lo = self.anneal(f_thresh, p_init, idx_lo, xgrid, p_min,
                                   **kwargs)

                self._vprint('Bounding above, limit = {}'.format(p_max))
                lmfit.confidence.restore_vals(orig, res.params)  # Reset
                for opstr, op in others_hi.items():  # Load other params
                    res.params[opstr].value = op
                p_hi = self.anneal(f_thresh, p_init, idx_hi, xgrid, p_max,
                                   **kwargs)
            else:
                self._vprint('Using naive grid bounds (no annealing!)')
                p_lo = xgrid[idx_lo]
                p_hi = xgrid[idx_hi]
        else:
            raise Exception('ERROR: fit type {} must be '
                            '\'full\' or \'simp\''.format(appr))

        # Restore values/chisqr in lmfit.Minimizer object
        lmfit.confidence.restore_vals(orig, res.params)
        res.params[pstr].vary = True
        res.chisqr = best_chisqr

        return p_lo, p_hi

    def prep_grid_anneal(self, chisqr_thresh, res, pstr):
        """Bracketing indices of chisqr_thresh in self.tab for pstr mesh

        This method only allows pstr = eta2, B0; i.e., application specific
        Output:
            xgrid for pstr, idx_lo/hi, pars_lo/hi
        """
        mu = res.params['mu'].value
        eta2_grid, B0_grid, fwhms_grid, chisqr_grid = self.grid_scan(mu)

        if pstr == 'eta2':
            xgrid = eta2_grid  # Already sorted
        elif pstr == 'B0':
            msk = np.argsort(B0_grid)
            xgrid = B0_grid[msk]
            eta2_grid = eta2_grid[msk]
            chisqr_grid = chisqr_grid[msk]
        else:
            raise Exception('Bad parameter name for full model errors')

        idx = np.searchsorted(xgrid, res.params[pstr].value)

        # Find chisqr crossings in pre-tabulated grid (may be incorrect)
        sgns = np.sign(chisqr_thresh(chisqr_grid))
        crossings = np.where(np.diff(sgns))[0]  # Left indices
        idx_lo = crossings[0] if any(crossings < idx) else 0
        # crossings[-1]+1 gets right index of crossing
        idx_hi = crossings[-1]+1 if any(crossings+1 > idx) else len(chisqr_grid)-1

        # Return values of "other" parameters @ lo/hi indices to ease fitting
        if pstr == 'eta2':
            pars_lo = {'B0': B0_grid[idx_lo]}
            pars_hi = {'B0': B0_grid[idx_hi]}
        elif pstr == 'B0':
            pars_lo = {'eta2': eta2_grid[idx_lo]}
            pars_hi = {'eta2': eta2_grid[idx_hi]}

        return xgrid, idx_lo, idx_hi, pars_lo, pars_hi

    def anneal(self, f, x0, idx, xgrid, x_lim, **kwargs):
        """Find root of f; search on/past grid from xgrid[idx_init] to x_lim

        Input:
            f is threshold function
            x0 is the best fit parameter value (not necessarily in xgrid)
            xgrid[idx] is where to start searching
            then searching outside if needed (up to x_lim)
            **kwargs passed to one_dir_root(...) (eps, verbose already set)
        Output:
            parameter bound for x0 (a root of f)
        """
        sgn = int(np.sign(x_lim - x0))
        self._vprint('Seeking error in {} direction;'.format(sgn),
                     'start: {:0.4e} (idx={}),'.format(xgrid[idx], idx),
                     'limit: {:0.4e}'.format(x_lim))
        if sgn == 0:
            self._vprint('Warning: sgn == 0 (?!), returning limit')
            return x_lim

        def get_eps(sgn_srch):
            """Safe forward/backwards step depending on sgn_src"""
            try:
                xstep = xgrid[idx + sgn_srch]
            except IndexError:
                xstep = xgrid[idx - sgn_srch]
            return abs(xgrid[idx] - xstep)

        # KEY ASSUMPTION:
        # if sgn < 0, then x_lim < xgrid[idx] < x0
        # if sgn > 0, then x0 < xgrid[idx] < x_lim

        # In practice, don't worry about x_lim
        if (sgn < 0 and xgrid[idx] < x0) or (sgn > 0 and x0 < xgrid[idx]):
            x_init = xgrid[idx]  # xgrid[idx] gives a better start position
        else:
            self._vprint('Searching outside grid')
            x_init = x0  # idx at grid edge

        f_init = f(x_init)

        if f_init < 0:
            self._vprint('Not crossed error bound, moving forward')
            return one_dir_root(f, x_init, x_lim, eps=get_eps(sgn),
                                **kwargs)
        elif f_init > 0:
            self._vprint('Found error crossing on grid, moving back')
            return one_dir_root(f, x_init, x0, eps=get_eps(-1*sgn),
                                **kwargs)
        else:
            self._vprint('A miracle has occurred, we are at crossing (?!)')
            return xgrid[idx]

    # -----------------------
    # Grid/table manipulation
    # -----------------------

    def grid_best(self, mu):
        """Find best fit to data/eps in (eta2, B0) grid of pretabulated FWHMs

        Input:
            data, eps (np.ndarray): data and errors to fit
            inds (np.ndarray): use only a subset of energy bands if desired
        Output:
            eta2 (float), B0 (float), fwhms (np.ndarray), chisqr (float)
        """
        eta2s, B0s, fwhms, chisqrs = self.grid_scan(mu)
        ind = np.argmin(chisqrs)
        return eta2s[ind], B0s[ind], fwhms[ind], chisqrs[ind]

    def grid_scan(self, mu):
        """Find best-fit B0s, FWHMs, chisqrs for eta2 grid

        Compute chisqr over (eta2, B0) grid for provided data/eps,
        subset with inds if necessary.  Cull best B0 value for each eta2.

        Input: data, eps (np.ndarray); eta2_dict (dict) grid of eta2, B0, FWHMs
        Output: SORTED eta2 values with best B0, fwhms, chisqr
        """
        if mu not in self.tab:
            mu_dist = np.abs(np.array(self.tab.keys()) - mu)
            mu2 = self.tab.keys()[mu_dist.argmin()]
            self._vprint('WARNING: mu ({:0.3f}) not in table; '
                         'using closest table value ({:0.3f}) '.format(mu,mu2))
            mu = mu2

        eta2_vals = np.sort(self.tab[mu].keys())
        B0_vals = np.array([])
        fwhm_vals = []  # Easier to build up, but not efficient
        chisqr_vals = np.array([])  # note: ugly code, oh well

        for eta2 in eta2_vals:
            B0_pts, fwhms = self.tab[mu][eta2]
            chisqr_pts = self.fwhm_scan(fwhms)

            # Previously tried interpolating better B0 fit, now scrapped
            ind_min = np.argmin(chisqr_pts)
            B0_vals = np.append(B0_vals, B0_pts[ind_min])
            fwhm_vals.append(fwhms[ind_min])
            chisqr_vals = np.append(chisqr_vals, chisqr_pts[ind_min])

        return eta2_vals, B0_vals, np.array(fwhm_vals), chisqr_vals

    def fwhm_scan(self, fwhms):
        """Compute chisqr values over pre-computed FWHMs (from B0 mesh)
        Uses stored data, eps, inds automatically
        """
        if self.inds is not None:  # Subset if needed
            fwhms = fwhms[:, self.inds]
        return map(lambda x: models.chi_squared(self.data, self.eps, x), fwhms)

    # ------------------------
    # Parameter space checking
    # ------------------------

    def check_eta2_grid(self, mus, fmts):
        """Plot chisqr-eta2 space from (eta2,B0) grids for multiple mu values
        JACKSON POLLOCK DOES ASTROPHYSICS

        Input:
            mus: list of mu values to check/plot
            fmts: list of linespec strings for plots
        Output:
            matplotlib.Axis object
        """
        ax = plt.gca()
        for mu, fmt in zip(mus, fmts):
            eta2s, B0s, fwhms, chisqrs = self.grid_scan(mu)
            ax.plot(eta2s, chisqrs, fmt)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(tuple(r'$\mu={:0.2f}$'.format(mu) for mu in mus), loc='best')
        fplot(r'$\eta_2$', r'$\chi^2$', axargs='tight')
        return ax

    def check_B0_grid(self, eta2, mus, fmts):
        """Plot chisqr(B0) from grid for fixed eta2, multiple mu values

        Verify that we can resolve chisq minimum in B0 grid

        Input:
            eta2: where in grid shall we slice?
            mus: list of mu values to check/plot
            fmts: list of linespec strings for plots
        Output:
            matplotlib.Axis object
        """
        ax = plt.gca()
        for mu, fmt in zip(mus, fmts):
            B0_vals, fwhms = self.tab[mu][eta2]
            chisqr_vals = self.fwhm_scan(fwhms)
            ax.plot(B0_vals*1e6, chisqr_vals, fmt)
        fplot(r'$B_0$ ($\mu$G)', r'$\chi^2$', axargs='tight')
        return ax

    def check_eta2_simp(self, eta2_vals, mus, fmts):
        """Plot chi-squared space for eta2, based on simple model"""
        ax = plt.gca()

        for mu, fmt in zip(mus, fmts):
            b0_vals = []
            chisqr_vals = []
            for eta2 in eta2_vals:
                res = self.fitter_simp(mu, eta2=eta2, eta2_free=False)
                b0_vals.append(res.params['B0'])
                chisqr_vals.append(res.chisqr)

            ax.plot(eta2_vals, chisqr_vals, fmt, lw=2)

        fplot(r'$\eta_2$ (-)', r'Best fit $\chi^2$', scales=('log','log'))
        ax.legend( tuple(r'$\mu = {:0.2f}$'.format(mu) for mu in mus),
                    loc='best')
        return ax


# ============================
# Confidence intervals, errors
# ============================
# Similar to `lmfit.conf_interval(...)`, but different metric for conf intvs

# NOTE could put inside Fitter? but, only instance variable used is verbose
# otherwise these are (almost) purely functional

# -----------------------------
# Formatting / wrangling output
# -----------------------------

def get_ci_errors(ci_dict, pstr, ci_val=0.683):
    """Read out +/- errors from dict of CI bounds (output in +/- order)"""
    x_bnds = ci_dict[pstr]
    x_lo, x_hi = np.sort([tup[1] for tup in x_bnds if tup[0] == ci_val])
    x = [tup[1] for tup in x_bnds if tup[0] == 0.][0]
    return x_hi - x, x - x_lo

def build_ci_dict(res, f_bounds, ci_vals=(0.683, 0.9)):
    """Confidence limits for free parameters in res (lmfit.Minimizer)

    Input:
        res is the output from lmfit.minimize(...), ALREADY MINIMIZED
        f_bounds is a function that takes CI, res, pstr; returns low/hi bounds
    Output:
        dict of ci bounds, structured like lmfit.conf_interval(...) output
    """
    ci2 = {}

    for pstr in [p for p in res.params if res.params[p].vary]:
        ci2[pstr] = [(0., res.params[pstr].value)]

        for ci_val in ci_vals:
            p_lo, p_hi = f_bounds(ci_val, res, pstr)
            ci2[pstr].append((ci_val, p_lo))
            ci2[pstr].append((ci_val, p_hi))

        ci2[pstr].sort(key=lambda x: x[1])

    return ci2

# --------------
# Finding bounds
# --------------

def stubborn_brentq(f, a2, a, b, b2, **kwargs):
    """brentq that tries 2x; if 2nd attempt fails you deal with it"""
    try:
        x = sp.optimize.brentq(f, a, b, **kwargs)
    except ValueError:
        print ('Inconsistent fit behavior; previously good values no '
               'longer bracket threshold. Making one more attempt...')
        x = sp.optimize.brentq(f, a2, b2, **kwargs)
    return x

def one_dir_root(f, x_init, x_lim, eps=None, adapt=True, verbose=False,
                 **kwargs):
    """Find root in function f from starting x with max/min limit
    Naive, brute force search in one direction, set by sgn(x_lim - x_init)
    If zero crossing is found, uses scipy.optimize.brentq to get root

    Designed specifically for chisqr bounds to compute confidence intervals
    Assumes monotonicity: if f(x_lim) doesn't cross, then you're SOL.
    But, code intended only for simple chisqr cases.  Beware of using elsewhere

    Adaptive search: if more than +/-1 from threshold, take 1.5x step size if
            % change in f (chisqr) < 50%, take /1.5 if % change > 75%
        If within +/-1 of threshold, same as above except use /absolute/ change

        Want to pass threshold and let brentq (with rtol specified) take over

    Input: f has call signature f(x)
           extra kwargs get sent to "stubborn brentq"
    Output: root within 2x machine precision, else x_lim if root was not found
    """
    if eps is None:  # Default is 1% of gap between initial/limit x
        eps = abs(0.01 * (x_lim - x_init)) * np.sign(x_lim - x_init)
    else:
        if abs(eps) > abs(x_lim - x_init):  # eps too big will overshoot!
            eps = abs(x_lim - x_init) * 0.25  # at least 4 steps to limit
        eps = abs(eps) * np.sign(x_lim - x_init)

    x = x_init
    dist = f(x_init)
    start_sgn = np.sign(dist)  # Allow crossings in either direction
    if start_sgn == 0:  # Edge case: we started right on crossing
        return x_init

    # Iterate until crossing is found (sign change)
    while np.sign(dist) == start_sgn:
        if x == x_lim:
            break
        elif verbose:
            print '\tIncrementing x, step size {}'.format(eps)

        prev_x = x
        prev_dist = dist

        x += eps
        if (eps > 0 and x > x_lim) or (eps < 0 and x < x_lim):
            x = x_lim
        dist = f(x)

        # Update step size
        # Note: since we can set brentq rtol, goal is to pass threshold
        # asap and let brentq take over (no need for full fits at too many
        # points.  So we can be somewhat aggressive in stepping.
        if adapt:
            print '\t  Prev chisqr dist {}, current {}'.format(prev_dist, dist)

            if dist > 1:  # Get close to threshold first
                if abs(dist - prev_dist)/dist < 0.5:
                    eps *= 1.5
                elif abs(dist - prev_dist)/dist > 0.75:
                    eps /= 1.5
            else:  # We want to pass threshold -- if not too far
                if abs(dist - prev_dist) < 0.5:
                    eps *= 1.5
                elif abs(dist - prev_dist) > 0.75:
                    eps /= 1.5

    # Broke out of loop (x = x_lim, dist < 0),
    # or hit crossing right on (dist == 0)
    if np.sign(dist) == start_sgn:  # Didn't pass crossing
        if verbose:
            print '\tCrossing not found, hit limit at {}'.format(x)
        return x
    elif verbose:
        print '\tCrossing found, attempting to bracket'

    # Brentq limits (first try [prev_x, x], then [x_lo, x_hi])
    x_lo = prev_x - 2*eps
    if (eps > 0 and x_lo < x_init) or (eps < 0 and x_lo > x_init):
        x_lo = x_init
    x_hi = x + eps

    lims = [x_lo, prev_x, x, x_hi]
    if eps < 0:
        lims.reverse()

    try:
        x_ret = stubborn_brentq(f, *lims, **kwargs)
    except ValueError:
        print ('WARNING: brentq failed; reporting guess; '
               'limit = {}, last step size eps={}'.format(x, eps))
        x_ret = x

    return x_ret


if __name__ == '__main__':
    main()

