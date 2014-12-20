"""
Functions to call models, fit models to data, and tabulate model FWHMs

Aaron Tran
2014 July 21
"""

from __future__ import division

import cPickle as pickle
from datetime import datetime
import inspect
import numpy as np
from operator import add as op_add
from operator import mul as op_mul
import os
import sys

import lmfit

import FullEfflength_port as fmp  # Full model code, Python port
from snr_catalog import SYNCH_B, SYNCH_CM, CD, NUKEV, BETA  # Constants

def main():
    pass

# ===========
# Chi squared
# ===========

# Kinda silly but multiple methods call this already.

def chi_squared(y, y_err, y_model):
    """Compute chi-squared statistic"""
    return np.sum( ((y_model - y) / y_err)**2 )

# =========================================
# Tabulate FWHM values from full model code
# =========================================

# If I ever change table format... make an n-d numpy array
# w/ coordinate axes (think meshgrid style)
# coords = {'eta2':[0, 1, 2, ... 1e5], 'B0':[100e-6, ..., 1e-3], ...}
# tab = np.ndarray([[[[....],[....],[...]]]])
# Then I can freely slice and dice.  Put NaNs wherever values aren't sampled.
# But, not necessary now.

# This could be split into a separate module, and make use of Fitter() and
# various convenience methods/wrappers.  Low priority.

def merge_tables(*args):
    """NOT IMPLEMENTED.  Lower priority at this time,
    but would be a nice feature to have.
    """
    # BEWARE NEED TO DEAL WITH REPEATED ETA2 VALUES...
    # ... well that's not so hard just merge the values and look for conflicts
    # in B0.  Only in this split case of Tycho will you need to address this
    # explicitly, and choose one or the other to discard.  Keep the one with
    # smaller rminarc.  Could make this interactive and prompt the user, or
    # automatically opt for the data w/ smaller rminarc.
    # BUT, rminarc is not stored as metadata with the tables.  ARGH.
    #.........
    # deal with this later.
    pass

def maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0,
    fname=None, f_B0_init=1.1, f_B0_step=0.15, **kwargs):
    """Tabulate FWHM values for set of SNR parameters in parameter space of
    (mu, eta2, B0), given some energy bands.

    First, grid over eta2.
    For each eta2, get range of B0 values that give reasonable FWHMs.
    Save B0 values and corresponding model FWHMs.

    kwargs are passed to width_cont for calculations

    TODO: Verify that results are independent of input FWHMs (mins, maxs, etc)
        Generate two tables w/ slightly different FWHM inputs/ranges
        Fit procedure should give results consistent within error
    """
    # Start logging to files (search for errors...)
    if fname is None:
        fname = '{}_gen_{}_grid_{}-{}-{}.pkl'.format(snr.name,
                datetime.now().strftime('%Y-%m-%d'), len(mu_vals),
                len(eta2_vals), n_B0)

    fname_b = os.path.splitext(fname)[0]
    log = fname_b + '.log'
    errlog = fname_b + '.errlog'
    stdout = TeeStdout(log, 'w')
    stderr = TeeStderr(errlog, 'w')
    np.set_printoptions(precision=3)

    # Spew a bunch of configuration parameters
    print '\nTabulating full model code FWHMs for SNR: {}'.format(snr.name)
    print 'Started: {}'.format(datetime.now())

    print '\nSNR parameters are:'
    print snr.config_log()

    print '\nDefault full model (width_cont) settings are:'
    args, _, _, dfltargs = inspect.getargspec(width_cont)
    for v in zip(args[-len(dfltargs):], dfltargs):
        print '\t{} = {}'.format(*v)
    print "'None' usually indicates kwarg is obtained from SNR object"
    print 'Overriding the following full model settings with:'
    for v in kwargs.items():
        print '\t{} = {}'.format(*v)

    print '\nGridding parameters are:'
    print 'Num grid points in mu, eta2, B0: {}, {}, {}+'.format(len(mu_vals),
            len(eta2_vals), n_B0)
    print 'Mu values are {}'.format(mu_vals)
    print 'eta2 values are:'
    print eta2_vals
    print 'B0 limits (for initial fitting) are {}'.format(snr.par_lims['B0'])
    print 'Gridding with FWHM limits:'
    print 'min: {}'.format(data_min)
    print 'max: {}'.format(data_max)

    print '\nAdditional parameters:'
    print 'File output stem: {}'.format(fname)
    print 'f_B0_init = {}'.format(f_B0_init)
    print 'f_B0_step = {}'.format(f_B0_step) # This is getting ridiculous

    # Loop over mu, eta2 grid; build up mu_dict, eta2_dict
    mu_dict = {}
    for mu in mu_vals:
        mu_dict[mu] = eta2_dict = {}

        for eta2 in eta2_vals:
            headstr = '(mu, eta2) = ({:0.2f}, {:0.2f})'.format(mu, eta2)
            print '\n', headstr, '\n', '-' * len(headstr)

            # Use simple model fit to set initial guess for B0
            # This fails badly for eta2 = 0, so the value of B0 above
            # had better be decent as a fallback.
            dmid = (data_max + data_min)/2
            res = simple_fit(snr, kevs, dmid, np.ones(len(dmid)),
                             mu, eta2=eta2, B0=150e-6,
                             mu_free=False, eta2_free=False, B0_free=True)
            p = res.params

            # Do the heavy work of computing FWHMs for many B0 values
            eta2_dict[eta2] = maketab_gridB0(snr, p, kevs, data_min, data_max,
                                             n_B0,
                                             f_B0_init = f_B0_init,
                                             f_B0_step = f_B0_step,
                                             **kwargs)

            # Save data regularly, but keep tabulating if save fails...
            try:
                with open(fname, 'w') as fpkl:
                    pickle.dump(mu_dict, fpkl)
                print 'Wrote table to file.'
            except:
                print 'ERROR: could not save table to file'
                fname = fname + '-badfname-save'  # Ad hoc workaround...
                print 'On next pass, will try saving to {}'.format(fname)

    # Stop logger
    print 'Finished: {}'.format(datetime.now())
    del stdout
    del stderr

    return mu_dict

def maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot,
                   f_B0_init=1.1, f_B0_step=0.15, **kwargs):
    """Grid over B0 at fixed eta2, mu (passed in via pars)

    1. Check that initial guess for B0 gives FWHMs in range of fwhms_min/max
       If initial guess is bad, increase/decrease B0 and try again.
    2. Generate grid over many values of B0, starting from initial B0.
       Store values of B0 and FWHMs from complex model.

    Parameters object, pars, will be modified throughout this method.

    Inputs:
        snr (snrcat.SupernovaRemnant): contains remnant information
        pars (lmfit.Parameters): contains mu, eta2, initial B0 guess

        kevs (np.array): energy bands, same ones used to generate fwhms_min/max
        fwhms_min, fwhms_max: must be numpy arrays

        n_tot (float): minimum number of evenly spaced data points.  In
            practice, this sets the max interval between rescaled mean widths
        f_B0_init (float): factor to change B0_init, if bad init guess
            must be greater than 1
        f_B0_step (float): max span_intv(...) step as fraction of initial B0

        **kwargs: passed straight to width_cont(...)

    Outputs:
        list of B0 values, list of FWHM value lists (for each B0)
    """

    def rscale(x):  # Grid with rescaled width averaged over FWHMS, r in [0,1]
        return np.mean((x-fwhms_min)/(fwhms_max - fwhms_min))

    # Use f_rscale for gridding in lieu of width_cont
    def f_rscale(grid_B0):
        """Rescale model width function"""
        pars.add('B0', value=grid_B0)  # Vary B0, other parameters same
        fwhms = width_cont(pars, kevs, snr, **kwargs)
        print '\t\tModel fwhms = {}'.format(fwhms)
        return rscale(fwhms), fwhms

    # Get FWHMs from initial guess, reinitialize B0 if initial guess is bad
    # WARNING: may fail if function is ill-behaved or f_B0_init is too large
    print 'Checking initial guess for B0'
    r_init, fwhms_init = f_rscale(pars['B0'].value)
    B0_init = pars['B0'].value
    if all(fwhms_init > fwhms_max):
        pars.add('B0', value=B0_init * f_B0_init)
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot,
                              f_B0_init, f_B0_step, **kwargs)
    elif all(fwhms_init < fwhms_min):
        pars.add('B0', value=B0_init / f_B0_init)
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot,
                              f_B0_init, f_B0_step, **kwargs)
    print 'Initial guess accepted'

    # Grid evenly over rescaling of FWHMS (rscale)
    n_thin = int(np.around(r_init * n_tot))  # Just for logging/debugging
    n_wide = int(np.around((1 - r_init) * n_tot))  # ditto
    dr = 1.0 / n_tot  # y-coord for span_intv

    print 'Using initial B0 value {} muG'.format(B0_init*1e6)
    print 'Now finding B0 values with FWHMs in range.'
    print 'Require {} values above, {} below initial B0'.format(n_thin, n_wide)

    # Calculate B0 values.  Only want values within or near r=[0,1], roughly
    if r_init < 1:
        print 'Computing values below initial B0'  # Increasing r towards 1
        B0_wide, r_wide, fwhms_wide = span_intv(B0_init, fwhms_init, r_init, 1,
                                                dr, f_rscale,
                                                dx_max = B0_init * f_B0_step)
    else:
        B0_wide, r_wide, fwhms_wide = [B0_init], [r_init], [fwhms_init]
    if r_init > 0:
        print 'Computing values above initial B0'  # Decreasing r towards 0
        B0_thin, r_thin, fwhms_thin = span_intv(B0_init, fwhms_init, r_init, 0,
                                                dr, f_rscale,
                                                dx_max = B0_init * f_B0_step)
    else:
        B0_thin, r_thin, fwhms_thin = [B0_init], [r_init], [fwhms_init]
    # Each call to span_intv first computes local derivative at the same spot
    # So user sees one repeated function call, but it's not stored in output

    # Combine lists
    B0_all = np.append(B0_thin, B0_wide)
    fwhms_all = np.append(fwhms_thin, fwhms_wide, axis=0)
    r_all = np.append(r_thin, r_wide)
    # Remove doubly included B0_init/r_init and sort
    B0_all, uniqsrt = np.unique(B0_all, return_index=True)
    fwhms_all = fwhms_all[uniqsrt]
    r_all = r_all[uniqsrt]

    # Select only FWHMs without rminarc errors for the "filling in FWHM range"
    # step, so that it does not get stuck at the discontinuity in r-values
    msk = np.all(fwhms_all < snr.rsarc*0.99, axis=1)  # Collapse FWHMs
    B0_all = B0_all[msk]
    fwhms_all = fwhms_all[msk]
    r_all = r_all[msk]

    # Fill in the gaps
    print 'Filling in FWHM range to achieve desired spacing'
    B0_all, fwhms_all, r_all = mind_the_gaps(B0_all, fwhms_all, r_all, dr,
                                             f_rscale)
    print 'Computed FWHMS for {} values of B0'.format(len(B0_all))

    return B0_all, fwhms_all  # No longer need r-values

def span_intv(x0, yaux0, y0, y_final, dy_step, f, dx_max=float('inf'),
             epsfcn_init=1e-2):
    """Find pts (x,y) with y-values spanning interval [y0, y_final]
    where y, yaux = f(x)

    As written, yaux (auxiliary function output) is specifically a list of
    numeric values.  Other auxiliary outputs may not work.

    This is Newton-Raphson iteration, but keeping intermediate values.
    dy_step is just a guideline for stepping torwards y_final
    To enforce data spacing, use mind_the_gaps

    WARNINGS: function must be monotonic, well-behaved, and have no extrema
    (i.e., derivative strictly nonzero and single-signed).
    There is no upper bound on the number of points generated.
    Code developed specifically for Sean's model fitting.

    Inputs:
        x0, yaux0, y0: initial data values (float, list, float)
        y_final (float): final y-value to iteratively reach
        dy_step (float): step size to use in y-coordinate
        f: function with call signature, output f(x) = y, yaux
        dx_max (float): maximum step size (dx), prevents overshoot
        epsfcn_init (float): fractional step size for initial derivative
            (to start the iteration), with dx = epsfcn_init * x0
    Outputs:
        x_vals, y_vals, yaux_vals from iterative process, with y_vals
        spanning the interval [y0, y_final].  Lists are not guaranteed
        to be sorted / ordered, but they may be so if the function f is
        well behaved.
    """
    sgn = np.sign(y_final - y0)
    dy_step = sgn * abs(dy_step)  # Sign of dy_step, set by y0/y_final
    dx_max = abs(dx_max)  # Force dx_max > 0

    # Convenience function to estimate step size dx, to reach dy_step
    def next_step(jcbn):
        dx = dy_step/jcbn  # if dy_step < 0 and jcbn < 0, dx > 0...
        return dx if abs(dx) < dx_max else np.sign(dx) * dx_max

    # Initialize storage for loop
    x_vals = np.array([x0])
    y_vals = np.array([y0])
    yaux_vals = np.array([yaux0])
    x_prev, y_prev = x_curr, y_curr = x0, y0  # Just for starters

    while sgn * y_curr < sgn * y_final:  # ad hoc stopping threshold
        # if dy_step > 0, proceed if y_curr < y_final
        # if dy_step < 0, proceed if y_curr > y_final
        if len(x_vals) > 1:
            jcbn = (y_curr - y_prev) / (x_curr - x_prev)  # Approx first deriv
        else:  # Initialize first deriv. if only (x0, y0) known
            dx_init = epsfcn_init * abs(x0)
            jcbn = (f(x0 + dx_init)[0] - y0)/dx_init

        dx = next_step(jcbn)

        # Shift last iteration values back
        x_prev, y_prev = x_curr, y_curr

        # Add new data point
        x_curr = x_prev + dx
        y_curr, yaux_curr = f(x_curr)

        x_vals = np.append(x_vals, x_curr)
        y_vals = np.append(y_vals, y_curr)
        yaux_vals = np.append(yaux_vals, [yaux_curr], axis=0)

    return x_vals, y_vals, yaux_vals

def mind_the_gaps(x, yaux, y, dy_max, f):
    """Fill data until y-spacing is always less than dy_max
    Assumptions: x, y, yaux are all sorted; y, yaux = f(x)

    Whenever a gap (dy > dy_max) is found, iteratively sample f at x-coordinate
    midpoints until all gaps are smaller than dy_max

    Input:
        x, y (np.array): data values to fill in
        yaux (np.array): auxiliary output of f
        dy_max (float): maximum y spacing, must be greater than zero
        f (function): function with call signature f(x).  Expect y, yaux = f(x)
    Output:
        x, y, yaux (np.array) with new data values
    """
    if dy_max <= 0:  # Prevent loop from running forever...
        raise ValueError('dy_max must be positive (input: {})'.format(dy_max))
    dy = np.abs(np.diff(y))
    while any(dy > dy_max):
        mask = dy > dy_max
        indgaps = np.where(mask)[0]

        x_midpts = (x[1:][mask] + x[:-1][mask]) / 2

        # Parse main + auxiliary output
        ytup_out = [f(xm) for xm in x_midpts]
        y_midpts = np.array([tup[0] for tup in ytup_out])
        yaux_midpts = np.array([tup[1] for tup in ytup_out])

        x = np.insert(x, indgaps + 1, x_midpts)
        y = np.insert(y, indgaps + 1, y_midpts)
        yaux = np.insert(yaux, indgaps + 1, yaux_midpts, axis=0)
        # Use axis=0 for lists as auxiliary output - namely lists of FWHMs
        dy = np.abs(np.diff(y))

    return x, yaux, y  # Match function call order...


# ======================================================
# Code to fit extracted image profiles to model profiles
# ======================================================

def profile_fit(snr, kevs, prfs, epss, r_prfs, mu, eta2=None, B0=None,
    r_trans=0, amp=1, multishift=False, ab_fit=None, mu_free=False, eta2_free=False,
    B0_free=True, rminarc_f = 1.2, model_kws=None, **lmfit_kws):
    """Fit measured thin rim profiles to model profiles, with single
    translation and amplitude shift for profiles at all energies
    (yet another convenience wrapper...)

    The model_kws 'get_prfs', 'irad_adapt', 'rminarc' are set/overriden so that
    profile fitting will work correctly.

    Inputs:
        snr, kevs as usual
        prfs, epss, r_prfs are lists of profiles, errors, r_grids
        mu, B0, eta2 (float) initial guesses
        mu_free, eta2_free (bool) which parameters shall vary in fits?
        ab_fit (float) if specified, allow ab to be a fit parameter
            (if you don't want it to vary, specify via model_kws instead)
            (BUT, enable idamp=True for this to work)
        **lmfit_kws takes extra kwargs for lmfit.minimize (includes lsq kws)

        model_kws={'icut':True, 'verbose':True, 'rminarc':60, ...}
    Output:
        lmfit.Minimizer with fit information / parameters
    """
    assert len(kevs) == len(prfs) == len(epss) == len(r_prfs)

    if model_kws is None:
        model_kws = {}
    model_kws['get_prfs'] = True
    model_kws['irad_adapt'] = False
    model_kws['get_fwhms'] = False
    model_kws['get_data'] = False

    # require rminarc > r_prfs at all energies, default safety factor 1.2
    rminarc = []
    for r_prf in r_prfs:
        rminarc.append(max(r_prf) - min(r_prf))
    model_kws['rminarc'] = rminarc_f * np.array(rminarc)

    p = lmfit.Parameters()
    for pstr in ['mu', 'B0', 'eta2']:
        pval = locals()[pstr]
        if pval is None:
            pval = snr.par_init[pstr]
        p.add(pstr, value=pval, vary=locals()[pstr+'_free'],
              min=snr.par_lims[pstr][0], max=snr.par_lims[pstr][1])

    if ab_fit is not None:
        p.add('ab', value=ab_fit, vary=True, min=0, max=1e3)

    if multishift:
        assert len(kevs) == len(r_trans) == len(amp)
        for n in xrange(len(amp)):
            p.add('r_trans_{:d}'.format(n), value=r_trans[n], vary=True)
            p.add('amp_{:d}'.format(n), value=amp[n], vary=True)
    else:
        p.add('r_trans', value=r_trans, vary=True)
        p.add('amp', value=amp, min=0, vary=True)

    res = lmfit.minimize(prf_objective_func, p,
                         args=(prfs, epss, r_prfs, multishift, kevs, snr),
                         kws=model_kws, **lmfit_kws)
    return res


def prf_objective_func(pars, prfs, epss, r_prfs, multishift, *args, **kwargs):
    """Objective function for profile fitting, for use with lmfit

    NOTES, please read:

    1. Code will not issue warning if best fit runs out of rminarc domain
       (interpolation of model profile just takes edge points in such cases).
       Please check that best fits lie within rminarc.

    2. for profile fitting to work correctly, you must supply the correct
       **kwargs (require profiles, don't get FWHMs, don't use adaptive rminarc,
       don't get useful "data". adaptive rminarc might work, but others are
       enforced)

    Inputs:
        pars: lmfit.Parameters object with entries B0, eta2, mu, ab,
              r_trans/amp or r_trans_[n]/amp_[n] (0 to len(prfs)-1)
              ab is optional depending on *args/**kwargs to width_cont
              r_trans/amp depend on multishift (see below)
        prfs: list of measured profiles (np.ndarray) corresponding to energies
              passed into kevs (*args; see width_cont docstring)
        epss: list of profile errors (np.ndarray), matching prfs/r_prfs
        r_prfs: list of r-coordinate grids (np.ndarray), matching prfs/epss
        multishift: indicate whether you want a single scaling/translation for
                    all profiles, or individual adjustments for each profile
                    pars entries must reflect your choice
        *args: arguments for width_cont (kevs, snr as of 2014 Nov 19)
        **kwargs: keyword arguments for width_cont, slew of flags goes here
    Output:
        np.ndarray of residuals for lmfit.minimize
    """
    assert len(prfs) == len(epss)
    assert len(epss) == len(r_prfs)

    if multishift:
        r_trans = [pars['r_trans_{:d}'.format(n)].value for n in xrange(len(prfs))]
        amp = [pars['amp_{:d}'.format(n)].value for n in xrange(len(prfs))]
    else:
        r_trans = [pars['r_trans'].value] * len(prfs)
        amp = [pars['amp'].value] * len(prfs)

    if 'verbose' not in kwargs or kwargs['verbose']:
        # Matches default of verbose=True in width_cont
        print '\tshifts (\") = {}'.format(r_trans)
        print '\tamps = {}'.format(amp)

    assert not kwargs['get_fwhms']
    assert 'get_data' not in kwargs or not kwargs['get_data']
    assert kwargs['get_prfs']

    r_prfs = map(op_add, r_prfs, r_trans)
    prfs = map(op_mul, prfs, amp)
    epss = map(op_mul, epss, amp)

    # Compute model profiles.  Don't need FWHMS (_)
    igrids, rgrids = width_cont(pars, *args, **kwargs)

    # Combine all profiles' weighted residuals
    res_all = []
    for _ in zip(prfs, epss, r_prfs, igrids, rgrids):
        prf, eps, r_prf, igrid, rgrid = _
        prf_fit = np.interp(r_prf, rgrid, igrid)
        res = (prf_fit - prf)/eps
        res_all.extend(list(res))

    return np.array(res_all)


# ==============================================
# Wrappers for full/simple model functions, fits
# ==============================================

# Input: snr, kevs, data/eps/mu/(eta2,B0);   output: best fit (mu)/eta2/B0

def simple_fit(snr, kevs, data, eps, mu, eta2=None, B0=None,
    mu_free=False, eta2_free=True, B0_free=True, **lmf_kws):
    """Perform a simple model fit (equation 6; Table 7 of Ressler et al.)
    A convenience wrapper for lmfit.minimize(objectify(width_dump), ...)

    Default is to fit both eta2, B0 at fixed mu.

    Inputs:
        kevs, data, eps (np.array) as usual
        mu, B0, eta2 (float) initial guesses, but mu fixed
        mu_free, eta2_free (bool) which parameters shall vary in fits?
        **lmf_kws takes extra kwargs for lmfit.minimize (includes lsq kws)
    Output:
        lmfit.Minimizer with fit information / parameters
    """
    p = lmfit.Parameters()
    for pstr in ['mu', 'B0', 'eta2']:
        pval = locals()[pstr]
        if pval is None:
            pval = snr.par_init[pstr]
        p.add(pstr, value = pval, vary = locals()[pstr+'_free'],
              min = snr.par_lims[pstr][0], max = snr.par_lims[pstr][1])

    res = lmfit.minimize(objectify(width_dump), p, args=(data, eps, kevs, snr),
                         **lmf_kws)
    return res

def full_fit(snr, kevs, data, eps, mu, eta2=None, B0=None, mu_free=False,
    eta2_free=True, B0_free=True, ab_fit=None, model_kws=None, **lmf_kws):
    """Perform full model fit (equation 12; Table 8 of Ressler et al.)
    Default fit: B0 free; mu, eta2 fixed

    Inputs:
        kevs, data, eps (np.ndarray) as usual
        mu, B0, eta2 (float) initial guesses, but mu fixed
        mu_free, eta2_free (bool) which parameters shall vary in fits?
        ab_fit (float) if specified, allow ab to be a fit parameter
            (if you don't want it to vary, specify via model_kws instead)
        **lmf_kws takes extra kwargs for lmfit.minimize (includes lsq kws)
        Probably should set epsfcn via **lsq_kws

        model_kws={'icut':True, 'verbose':True, 'rminarc':60, ...}
    Output:
        lmfit.Minimizer with fit information / parameters
    """
    if model_kws is None:
        model_kws = {}

    p = lmfit.Parameters()
    for pstr in ['mu', 'B0', 'eta2']:
        pval = locals()[pstr]
        if pval is None:
            pval = snr.par_init[pstr]
        p.add(pstr, value=pval, vary=locals()[pstr+'_free'],
              min=snr.par_lims[pstr][0], max=snr.par_lims[pstr][1])

    if ab_fit is not None:
        p.add('ab', value=ab_fit, vary=True, min=0, max=1e3)

    res = lmfit.minimize(objectify(width_cont), p, args=(data, eps, kevs, snr),
                         kws=model_kws, **lmf_kws)
    return res

# Input: snr, kevs, mu/eta2/B0; output: model FWHMs

def full_width(snr, kevs, mu, eta2, B0, **kwargs):
    """Full model FWHMs for given SNR, input parameters, energy bands
    (wrapper for width_cont, can pass through intensities/r-coord grids too)
    Input:
        snr (snrcat.SupernovaRemnant): container w/various physical parameters
        kevs (np.ndarray): energies at which to compute profile FWHMs
        mu, eta2, B0 (scalars): diffusion parameters, magnetic field
        **kwargs: passed to width_cont(...), many important twiddle-ables
    Output:
        np.ndarray of simple model FWHMs
    """
    p = lmfit.Parameters()
    p.add('mu', value=mu)
    p.add('eta2',value=eta2)
    p.add('B0', value=B0)
    return width_cont(p, kevs, snr, **kwargs)

def simple_width(snr, kevs, mu, eta2, B0):
    """Simple model FWHMs for given SNR, input parameters, energy bands
    Input:
        snr (snrcat.SupernovaRemnant): container w/various physical parameters
        kevs (scalar, np.ndarray): energies at which to compute profile FWHMs
        mu, eta2, B0 (scalars): diffusion parameters, magnetic field
    Output:
        np.ndarray of simple model FWHMs
    """
    p = lmfit.Parameters()
    p.add('mu', value=mu)
    p.add('eta2',value=eta2)
    p.add('B0', value=B0)
    return width_dump(p, kevs, snr)


# =============================
# Model functions, objectify(f)
# =============================
# Input: params, kevs, snr, **kwargs; output: model FWHMs

# TODO move various defaults into snr_catalog
def width_cont(params, kevs, snr, verbose=True, rminarc=None, icut=None,
    irmax=None, iradmax=None, ixmax=None, irad_adapt=True, irad_adapt_f=1.2,
    idamp=False, damp_ab=0.05, damp_bmin=5.0e-6, fgfname='fglists_mod.dat',
    itmax=200, inmax=50, irhomax=2000, get_prfs=False, get_data=False,
    get_fwhms=True):
    """Width function, wrapper for Python port of full model (equation 12)

    All **kwargs not specified (excepting verbose) are taken from snr object

    Inputs:
        params: lmfit.Parameters object with entries B0, eta2, mu
            B0 (Gauss), eta2 (-), mu (-)
            IF ab (-) is in params, it will override damp_ab (!)
            BUT, you still have to toggle idamp=True or it will do nothing
        kevs: scalar/vector of observed photon energies
        snr: object with SNR physical attributes
    **kwargs:
        rminarc: distance downstream of shock out to which to compute profiles
                 units arcsec, profile computed on [rsarc - rminarc, rsarc]
                 distance should be slightly greater than computed FWHMs
                 (and rsarc > rminarc, of course)
                 np.ndarray, int, or float
        icut: boolean; toggle energy cutoff in injected e- energy spectrum
        irmax: resolution on intensity profile (radial coordinate)
        iradmax: resolution on tabulated e- distribution (radial coordinate)
                 (energy resolution for e- distr. set by fglists.dat)
        ixmax: resolution on line of sight for emissivity integration

        irad_adapt: boolean, compute FWHMs 2x for better precision
                    (slower, but should be much more precise and is
                    comparatively faster than increasing grid resolution)
        irad_adapt_f: float, should be > 1; margin-of-error for rminarc in
                      adaptive calculation

        idamp: boolean; toggle magnetic damping
        damp_ab: float between [0,1], damping lengthscale
        damp_bmin: float, minimum magnetic field for damping

        fgfname: str, pre-compiled table of 1-particle synchrotron emissivity
                 originally from Pacholczyk

        itmax, inmax: internal integral resolutions for e- distributions
                      (Sean's defaults were itmax=1000, inmax=100)
        get_prfs: get intensity profiles and r grid values alongside FWHMs
                  NOTE -- this BREAKS the fitting functionality
                  so, obviously, don't use if you are fitting stuff
        get_data: get useful information about profile shape (min/max pos,
            existence of a trough, % above/below rim max)
        get_fwhms: True by default. Saves marginally on computation time if you
            don't want FWHMs!

    Outputs:
        np.array of modeled FWHMs (arcsec) for each energy in 'kevs'
        if get_prfs, then also returns intensity profiles and r grid values
        as 3-tuple of fwhms, intensities, rgrids
    """
    # ------------
    # Sanity check
    # ------------

    assert get_prfs or get_data or get_fwhms
    # You have to get SOMETHING out...

    # -----------------------
    # Load parameter/settings
    # -----------------------

    B0 = params['B0'].value
    eta2 = params['eta2'].value
    mu = params['mu'].value

    if 'ab' in params:
        damp_ab = params['ab'].value

    vs = snr.vs
    v0 = snr.v0()
    rs = snr.rs()
    rsarc = snr.rsarc
    s = snr.s

    if rminarc is None:  # Use default rminarc
        rminarc = snr.rminarc * np.ones(len(kevs))
    elif isinstance(rminarc, int) or isinstance(rminarc, float):
        rminarc = rminarc * np.ones(len(kevs))
    # else: supplied rminarc should be np.array

    if icut is None:  # Lulz...
        icut = snr.icut
    if irmax is None:
        irmax = snr.irmax
    if iradmax is None:
        iradmax = snr.iradmax
    if ixmax is None:
        ixmax = snr.ixmax

    # ------------------------
    # Begin actual computation
    # ------------------------

    if irad_adapt:
        iradmax_prelim = int(iradmax / 2.)  # Just a rough calculation

        out_prelim = fmp.fefflen(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s,
                                 rminarc, icut, irmax, iradmax_prelim, ixmax,
                                 idamp, damp_ab, damp_bmin, fgfname,
                                 itmax, inmax, irhomax,
                                 False, False,  # get_prfs / get_data
                                 True)  # get_fwhms
        fwhms_prelim = out_prelim[0]

        # Replace error-code values
        res_msk, box_msk = _check_calc_errs(fwhms_prelim, rminarc)
        nan_msk = np.logical_not(np.isfinite(fwhms_prelim))
        big_msk = np.logical_or(box_msk, nan_msk)
        # box resolution errors + any NaNs/infs (unlikely, but in case)
        if any(big_msk):
            fwhms_prelim[big_msk] = rminarc[big_msk]
        # Resolution errors, take the smallest fwhm
        if all(res_msk):  # Extreme edge case...
            fwhms_prelim = rminarc / (irad_adapt_f**2)  # Hacky solution
        elif any(res_msk):
            fwhms_prelim[res_msk] = min(fwhms_prelim[np.logical_not(res_msk)])

        # Construct adapted rminarc
        rminarc2 = irad_adapt_f * fwhms_prelim

    # Print settings, set rminarc to new adaptive value
    if verbose:
        if idamp:
            parstr = ('\tFull model: B0 = {:0.3f} muG; eta2 = {:0.3f}; '
                      'mu = {:0.3f}; ab = {:0.3f}, Bmin = {:0.3f} muG')
            parstr = parstr.format(B0*1e6, eta2, mu, damp_ab, damp_bmin*1e6)
            print parstr
        else:
            print ('\tFull model: B0 = {:0.3f} muG; eta2 = {:0.3f}; '
                   'mu = {:0.3f};').format(B0*1e6, eta2, mu)

        if irad_adapt:
            print '\tinit rminarc = {}; adapted = {}'.format(rminarc, rminarc2)
        else:
            print '\trminarc = {}'.format(rminarc)

    # Set rminarc to new adaptive value
    if irad_adapt:
        rminarc = rminarc2

    # Python full model port
    out = fmp.fefflen(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s,
                      rminarc, icut, irmax, iradmax, ixmax,
                      idamp, damp_ab, damp_bmin, fgfname,
                      itmax, inmax, irhomax, get_prfs, get_data, get_fwhms)

    if get_fwhms:
        fwhms = out[0]
        res_msk, box_msk = _check_calc_errs(fwhms, rminarc)
        if any(res_msk):
            print '\t\tResolution error! at {}'.format(kevs[res_msk])
        if any(box_msk):
            print '\t\tBox error! at {}'.format(kevs[box_msk])

    # Preserving old behavior for backwards compatability
    # BUT should be updated throughout codebase if possible (TODO)
    # to just accept out directly.
    if get_fwhms and not (get_prfs or get_data):
        return fwhms
    return out

def _check_calc_errs(fwhms, rminarc):
    """Check for resolution or box errors (as in Sean's code)"""
    reserr = fwhms < np.finfo(float).eps * 2  # Any fwhms <= 0 (error output)
    boxerr = fwhms > (rminarc - np.finfo(float).eps * 2)  # Any FWHMs = 1
    return reserr, boxerr

def width_dump(params, kevs, snr):
    """Width function from catastrophic dump model (equation 6), code taken
    from Sean's code Widthfun.py

    If eta2 is very small: l_ad / l_diff > 30, or D / (v0**2 * tsynch) < 1/900,
    use pure advective solution to avoid singularity (should ease fitting)
    Same criterion as in Sean's full model code

    Inputs:
        params: lmfit.Parameters object with entries B0, eta2, mu
            B0 (Gauss), eta2 (-), mu (-)
        kevs: scalar/vector of observed photon energies
        snr: object with SNR physical attributes
    Output:
        np.array of modeled FWHMs (arcsec) for each energy in 'kevs'
    """
    B0 = params['B0'].value
    eta2 = params['eta2'].value
    mu = params['mu'].value

    # SNR parameters
    v0 = snr.v0()
    rs = snr.rs()
    rsarc = snr.rsarc

    # "universal" constants
    b = SYNCH_B
    cm = SYNCH_CM
    Cd = CD
    nu2 = NUKEV * 2.0
    beta = BETA

    # eta = eta2 * (2 keV)**(1-mu)
    # eta in code maps to \eta*(E_h)^(1-\mu) in paper
    eta = eta2 * np.sqrt(nu2/(cm*B0))**(1-mu)
    nus = kevs * nu2/2.
    E = np.sqrt(nus/(cm*B0)) # Energy of electron radiating at freq. nu
                             # in magnetic field B0
    tsynch = 1.0/(b*E*B0**2)
    D = eta*Cd*E**mu/B0

    # Equation 6, width = beta * a
    # Note that D, tsynch, nus, E, etc.. are vectors
    width = beta * 2*D/v0 / (np.sqrt(1 + (4*D/v0**2/tsynch)) - 1)

    # If diffusion is small, use advection only solution
    # similar to full model code (here, prevents blowup in simple model fits)
    msk = D / (v0**2 / tsynch) < 1./900
    try:
        width[msk] = v0*tsynch[msk]
    except TypeError: # width is a scalar
        if msk:
            width = v0*tsynch

    return width*rsarc/rs  # Convert width (cm) to arcsec

def objectify(f):
    """Generate objective function g(params, data, eps_data, *args, **kwargs)
    from model function f(params, *args, **kwargs), for lmfit.minimize()

    params is a lmfit.Parameters() object
    *args typically takes the "x"-coordinates (abscissae) corresponding
    to the data you seek to fit
    **kwargs optional parameters, as usual

    Input:
        f: model function for data, signature f(pars, x)
    Output:
        objective function with additional arguments data, eps_data
            data: measured data to be fitted
            eps_data: measured data errors (used as weights)
    """
    return lambda pars, data, eps_data, *args, **kwargs: \
        (f(pars, *args, **kwargs) - data)/eps_data

# ===========================
# Prototype code (not usable)
# ===========================

#from __future__ import division
#import numpy as np
#import scipy as sp
#from scipy.integrate import odeint
#import matplotlib.pyplot as plt

#def width_dump_distr():
#    """Code to numerically solve for f(E,x) for simple model with magnetic
#    field damping.  See my notes from 2014 October 18.  Not functional at this
#    time (blows up for f'(x=0) != 0, basically)
#    If used, would be integrated into width_dump appropriately...
#    """
#
#    CD = 2.083e19   # c/(3e)
#    b = 1.57e-3     # Synchrotron cooling constant
#
#    eta = 1  # Regular eta, not eta2
#    mu = 1
#
#    v = 1.25e8  # shock veloc, 1/4 of 5e8 = 5000 km/s
#    E = 36.5    # e- energy for 1 keV photon at B = 100 \muG
#    Bmin = 5e-6
#    B0 = 30e-6
#    ab = 1e17   # Damping length ~1% of Tycho's shock radius
#
#    rs_tycho = 240 * np.pi/180. /3600. * 3.0 * 1e3 * 3.085678e16 * 1e2
#    print rs_tycho
#
#    def B(x):
#        """Damped B-field"""
#        return Bmin + (B0 - Bmin) * np.exp(-x/ab)
#
#    def Bprime(x):
#        """Spatial derivative of damped B-field"""
#        return -1/ab * (B0 - Bmin) * np.exp(-x/ab)
#
#    def D(x):
#        """Diffusion coefficient"""
#        return eta * CD * E**mu / B(x)
#
#    def odefunc(vec, x):
#        """Function for ODE derivatives at x """
#        f, g = vec
#        c1 = - b* np.power(B(x),2) * E / D(x)
#        c2 = v/D(x) + Bprime(x) / B(x)
#        return g, c1*f + c2*g
#
#    ic = [1, 0]  # f(x=0) = 1, f'(x=0) = 0
#
#    xvec = np.linspace(0,100*ab,1000)
#    svec = odeint(odefunc, ic, xvec); fvec = svec[:,0]
#
#    plt.plot(xvec,fvec)
#    plt.ylim(0,1)
#    plt.show()
#
#    return xvec, svec


# ===================================
# Useful tees for log files in Python
# ===================================

class TeeStdout(object):
    """From stackoverflow.com/a/616686 and Python mailing list
    mail.python.org/pipermail/python-list/2007-May/438106.html
    Not ideal (issues with __del__, close?)... but it works, so oh well.
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()

class TeeStderr(object):
    """Source as above. Too lazy to refactor code, and
    can't reassign sys.stderr when passed in anyways
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stderr = sys.stderr
        sys.stderr = self
    def __del__(self):
        sys.stderr = self.stderr
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stderr.write(data)
    def flush(self):
        self.file.flush()


if __name__ == '__main__':
    main()
