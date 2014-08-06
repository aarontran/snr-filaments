"""
Model fitting functions for measured FWHMs
Does NOT recompile f2py-generated fullmodel.so

Methods (for my own reference):

    # Fitting wrappers
    simple_fit(snr, kevs, data, eps, mu, eta2=1, B0=150e-6, ...)
    full_fitter(snr, kevs, data, eps, mu)
    full_fit(snr, kevs, data, eps, mu, eta2=1, B0=15e-6, ...)

    # Interactive fitting
    manual_fit(snr, kevs, data, eps)
        _get_float(prompt)

    # Tabulating
    merge_tables(...)
    maketab(...)
        maketab_gridB0(...)
        span_intv(...)
        mind_the_gaps(...)
        classes TeeStdout, TeeStderr

    # Functions to fit
    width_cont(...)
    width_dump(...)
    objectify(f)

    chi_squared(y, y_err, y_model)


Aaron Tran
2014 July 21
"""

from __future__ import division

import lmfit
import matplotlib.pyplot as plt
import numpy as np
from numpy import f2py
import scipy as sp
from scipy import optimize

import cPickle as pickle
from datetime import datetime
import os
import sys

import fullmodel as fm
fm.readfglists()  # Only do this once
from fplot import fplot
import snr_catalog as snrcat


def main():
    pass


# ==================================================
# Fit from just initial guesses (simple, full model)
# ==================================================

# Convenience wrappers for the fit functions,
# feed in parameters and get widths out quickly.

def simple_fit(snr, kevs, data, eps, mu, eta2=1, B0=150e-6,
               mu_free=False, eta2_free=True, B0_free=True):
    """Perform a simple model fit (equation 6; Table 7 of Ressler et al.)
    A convenience wrapper for lmfit.minimize(objectify(width_dump), ...)

    Default is to fit both eta2, B0 at fixed mu.
    
    Inputs:
        kevs, data, eps (np.array) as usual
        mu, B0, eta2 (float) initial guesses, but mu fixed
        mu_free, eta2_free (bool) which parameters shall vary in fits?
    Output:
        lmfit.Minimizer with fit information / parameters
    """
    p = lmfit.Parameters()
    p.add('mu', value=mu, vary=mu_free)
    p.add('B0', value=B0, min=1e-6, max=1e-2, vary=B0_free)
    p.add('eta2', value=eta2, min=5e-16, max=1e5, vary=eta2_free)
    # Must be nonzero to prevent width_dump from going singular
    res = lmfit.minimize(objectify(width_dump), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq')
    return res


def full_fitter(snr, kevs, data, eps, mu):
    """Wrapper function to package SNR and data dependence"""
    return lambda **kwargs: full_fit(snr, kevs, data, eps, mu, **kwargs)


def full_fit(snr, kevs, data, eps, mu, eta2=1, B0=150e-6,
             mu_free=False, eta2_free=False, B0_free=True,
             rminarc=None, verbose=True, **lsq_kws):
    """Perform a full model fit (equation 12; Table 8 of Ressler et al.)
    A convenience wrapper for lmfit.minimize(objectify(width_cont), ...)

    Default is to fit B0 at fixed mu, eta2.

    Inputs:
        kevs, data, eps (np.array) as usual
        mu, B0, eta2 (float) initial guesses, but mu fixed
        mu_free, eta2_free (bool) which parameters shall vary in fits?
        Probably should set epsfcn via **lsq_kws
    Output:
        lmfit.Minimizer with fit information / parameters
    """
    p = lmfit.Parameters()
    p.add('mu', value=mu, vary=mu_free)
    p.add('B0', value=B0, min=1e-6, max=1e-2, vary=B0_free)
    p.add('eta2', value=eta2, min=5e-16, max=1e5, vary=eta2_free)

    # TODO: allow kws to width_cont?... function to get kwargs dict
    res = lmfit.minimize(objectify(width_cont), p,
                         args=(data, eps, kevs, snr),
                         kws={'icut':1, 'verbose':verbose, 'rminarc':rminarc},
                         # width_cont settings
                         method='leastsq', **lsq_kws)
    return res


# ============================================
# Interactive manual fitting (full model code)
# ============================================

def manual_fit(snr, kevs, data, eps):
    """Interactively prompts user for values of B0, eta2, mu
    Prints FWHMs, residuals, and chi-squared values in the process
    And, plots the data + fitted model
    Also prints m_E, computed point to point (just for comparison)

    Note: doesn't really deal with bad input

    Input: snr object and data (x, y, eps) to be fitted
    Output: lmfit.Parameters object with best obtained fit values
    """

    print 'Manual fitting routine: enter q to print best fit and quit.\n'

    chisq_best = float('inf')
    w_model_best = None
    p_best = None

    while True:
        B0 = _get_float('Enter B0 (G): ')
        eta2 = _get_float('Enter eta2 (-): ')
        mu = _get_float('Enter mu (-): ')

        if any(np.isnan([B0, eta2, mu])):
            break  # Received exit/quit request

        p = lmfit.Parameters()
        p.add('B0', value=B0)
        p.add('eta2', value=eta2)
        p.add('mu', value=mu)

        adj = _get_float('Change model settings? Enter 0/1 for no/yes: ')
        if adj == 0:
            w_model = width_cont(p, kevs, snr)
        else:
            rminarc = _get_float('rminarc (default: 60): ')
            icut = _get_float('icut (default: 1): ')
            # Could add option to change resolution too
            w_model = width_cont(p, kevs, snr, rminarc=rminarc, icut=icut)

        chisq = chi_squared(data, eps, w_model)
        wresid = (w_model - data)/eps

        print ('\nObsvd fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*data)
        print ('Model fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*w_model)
        print ('Wghtd resid: ' + len(kevs)*'{:0.3f}, ')[:-2].format(*wresid)
        print 'Chi^2: {:0.3f}'.format(chisq)

        plt.clf()
        plt.errorbar(kevs, data, eps, fmt='bo')
        plt.plot(kevs, w_model, '-k.')
        plt.draw()
        plt.show(block=False)

        if chisq < chisq_best:
            print '\nImproved fit, saving parameters'
            chisq_best = chisq
            w_model_best = w_model
            p_best = p
        else:
            print '\nBest fit so far:'
            print ('Best fwhms:  ' +
                   len(kevs)*'{:0.2f}, ')[:-2].format(*w_model_best)
            print 'Chi^2: {:0.3f}'.format(chisq_best)

        print '--------------------------------'
    
    if p_best is not None:
        plt.close()
        print '\nDone with manual fit. Best fit parameters are:'
        for key in p_best.keys():
            print '{} = {:g}'.format(p_best[key].name, p_best[key].value)

    return p_best


def _get_float(prompt):
    """Prompt user to input float, checking for exits / bad floats"""
    while True:
        try:
            uinput = raw_input(prompt).strip()
            if uinput in ['q', 'quit', 'exit']:
                uinput = float('NaN') # Pass NaN to indicate exit/quit
            else:
                uinput = float(uinput)
            break
        except ValueError:
            print '\nInvalid float, try again (enter q to quit)\n'
    return uinput


def chi_squared(y, y_err, y_model):
    """Compute chi-squared statistic"""
    return np.sum( ((y_model - y) / y_err)**2 )


# =========================================
# Tabulate FWHM values from full model code
# =========================================

def merge_tables(*args):
    """Does what it says.  Lower priority at this time,
    but would be a nice feature to have."""
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
            fname=None, rminarc=None, f_rminarc=1.2, f_B0_init=1.1,
            f_B0_step=0.15):
    """Tabulate FWHM values for set of SNR parameters in parameter space of
    (mu, eta2, B0), given some energy bands.

    First, grid over eta2.
    For each eta2, get range of B0 values that give reasonable FWHMs.
    Save B0 values and corresponding model FWHMs.

    TODO: Verify that results are independent of input FWHMs (mins, maxs, etc)
        Generate two tables w/ slightly different FWHM inputs/ranges
        Fit procedure should give results consistent within error

    Think about full model fit errors?  Once we have best fit B0, eta2, mu for a
    given region/filament, then deal with error.  Here we're just making
    a table of possible/candidate values.
    """

    # Mildly convenient to define rminarc here
    if rminarc is None:
        rminarc = data_max * f_rminarc

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
    np.set_printoptions(precision=2)

    # Spew a bunch of configuration parameters
    print '\nTabulating full model code FWHMs for SNR: {}'.format(snr.name)
    print '\nStarted: {}'.format(datetime.now())
    print 'SNR parameters are (rminarc will be modified):'
    for v in vars(snr):
        print '\t{} = {}'.format(v, vars(snr)[v])
    print '\tDerived SNR parameters:'
    print '\tv0 = {}'.format(snr.v0())
    print '\trs = {}'.format(snr.rs())

    print '\nGridding parameters are:'
    print 'Resolution in mu, eta2, B0: {}, {}, {}+'.format(len(mu_vals),
            len(eta2_vals), n_B0)
    print 'Mu values are {}'.format(mu_vals)
    print 'eta2 values are {}'.format(eta2_vals)
    print 'Gridding with FWHM limits:'
    print 'min: {}'.format(data_min)
    print 'max: {}'.format(data_max)
    print 'rminarc: {}'.format(rminarc)

    print '\nAdditional parameters:'
    print 'File output stem: {}'.format(fname)
    print 'f_rminarc = {}'.format(f_rminarc)
    print 'f_B0_init = {}'.format(f_B0_init)
    print 'f_B0_step = {}'.format(f_B0_step) # This is getting ridiculous

    # Loop over mu, eta2 grid
    mu_dict = {}
    for mu in mu_vals:
        mu_dict[mu] = eta2_dict = {}
        p = lmfit.Parameters()
        p.add('mu', value=mu, vary=False)

        for eta2 in eta2_vals:
            print '\n(mu, eta2) = ({:0.2f}, {:0.2f})'.format(mu, eta2)
            print '---------------------------------'

            p.add('eta2', value=eta2, vary=False)
            p.add('B0', value=150e-6, min = 1e-6, max=1e-2)

            # Use simple model fit to set initial guess for B0
            # This does fail badly for eta2 = 0, so the value of B0 above
            # had better be decent as a fallback.
            dmid = (data_max + data_min)/2
            res = lmfit.minimize(objectify(width_dump), p, method='leastsq',
                                 args=(dmid, np.ones(len(dmid)), kevs, snr))

            # Do the heavy work of computing FWHMs for many B0 values
            eta2_dict[eta2] = maketab_gridB0(snr, p, kevs, data_min, data_max,
                                             rminarc, n_B0, f_B0_init,
                                             f_B0_step)

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


def maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, rminarc, n_tot,
                   f_B0_init=1.1, f_B0_step=0.15):
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
        fwhms_min, fwhms_max, rminarc: must be numpy arrays

        n_tot (float): minimum number of evenly spaced data points.  In
            practice, this sets the max interval between rescaled mean widths
        f_B0_init (float): factor to change B0_init, if bad init guess
            must be greater than 1
        f_B0_step (float): max span_intv(...) step as fraction of initial B0

    Outputs:
        list of B0 values, list of FWHM value lists (for each B0)
    """

    def rscale(x):  # Grid with rescaled width averaged over FWHMS, r in [0,1]
        return np.mean((x-fwhms_min)/(fwhms_max - fwhms_min))

    # Use f_rscale for gridding in lieu of width_cont
    def f_rscale(grid_B0):
        """Rescale model width function"""
        pars.add('B0', value=grid_B0)  # Vary B0, other parameters same

        fwhms = width_cont(pars, kevs, snr, rminarc = rminarc)

        if any(fwhms <= 1e-5):
            print "Resolution error!"
        if any(fwhms >= snr.rsarc - 1e-5):
            print "Box length error!"
        print '\tModel fwhms = {}'.format(fwhms)

        return rscale(fwhms), fwhms

    # Get FWHMs from initial guess, reinitialize B0 if initial guess is bad
    # WARNING: may fail if function is ill-behaved or f_B0_init is too large
    print 'Checking initial guess for B0'
    r_init, fwhms_init = f_rscale(pars['B0'].value)
    B0_init = pars['B0'].value
    if all(fwhms_init > fwhms_max):
        pars.add('B0', value=B0_init * f_B0_init)
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot,
                              f_B0_init, f_B0_step)
    elif all(fwhms_init < fwhms_min):
        pars.add('B0', value=B0_init / f_B0_init)
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot,
                              f_B0_init, f_B0_step)  # KEEP UPDATED... ugh
    print 'Initial guess accepted'

    # Grid evenly over rescaling of FWHMS (rscale)
    n_thin = int(np.around(r_init * n_tot))  # Just for logging/debugging
    n_wide = int(np.around((1 - r_init) * n_tot))  # ditto
    dr = 1.0 / n_tot  # y-coord for span_intv

    print 'Using initial B0 value {} muG'.format(B0_init*1e6)
    print 'Now finding B0 values with FWHMs in range.'
    print 'Require {} values above, {} below initial B0'.format(n_thin, n_wide)

    # Now, actually calculate B0 values.  Only want values in/near r=[0,1]
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

    # Fill in the gaps
    print 'Filling in FWHM range to achieve desired spacing'
    B0_all, fwhms_all, r_all = mind_the_gaps(B0_all, fwhms_all, r_all, dr,
                                             f_rscale)

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


# =============================
# Model functions, objectify(f)
# =============================

def width_cont(params, kevs, snr, verbose=True, rminarc=None, icut=1,
               irmax=400, iradmax=100, ixmax=500):
    """Width function, wraps numerical computation in fm.fullefflengthsub
    which solves continuous loss model (equation 12)

    Be careful -- rminarc depends on snr being fitted.
    If not provided, rminarc is set by an SNR "default"
    (e.g., SN 1006 has default rminarc=60)

    Inputs:
        params: lmfit.Parameters object with entries B0, eta2, mu
            B0 (Gauss), eta2 (-), mu (-)
        kevs: scalar/vector of observed photon energies
        snr: object with SNR physical attributes
    **kwargs:
        rminarc: distance downstream of shock out to which to compute profiles
                 units arcsec, profile computed on [rsarc - rminarc, rsarc]
                 distance should be slightly greater than computed FWHMs
                 (and rsarc > rminarc, of course)
                 either np.ndarray, or a numeric type
        icut: 0 or 1; toggle electron spectrum cutoff for f(E,x)
        irmax: resolution on intensity profile (radial coordinate)
        iradmax: resolution on tabulated e- distribution (radial coordinate)
                 (energy resolution for e- distr. set by fglists.dat)
        ixmax: resolution on line of sight for emissivity integration
    Outputs:
        np.array of modeled FWHMs (arcsec) for each energy in 'kevs'
    """
    B0 = params['B0'].value
    eta2 = params['eta2'].value
    mu = params['mu'].value

    vs = snr.vs
    v0 = snr.v0()
    rs = snr.rs()
    rsarc = snr.rsarc
    s = snr.s

    # Case handling for rminarc, keep old methods from breaking
    def is_num(x):  # SUPER SKETCH BUT OH WELL
        return isinstance(x, int) or isinstance(x, float)

    if rminarc is None:
        rminarc = snr.rminarc
        if is_num(snr.rminarc):
            rminarc = rminarc * np.ones(len(kevs))
    elif is_num(rminarc)
        rminarc = rminarc * np.ones(len(kevs))
    # else: rminarc is already an array

    if verbose:
        print ('\tFunction call with B0 = {:0.3f} muG; eta2 = {:0.3f}; '
               'mu = {:0.3f}; rminarc = {}').format(B0*1e6, eta2, mu, rminarc)

    # FORTRAN subroutine call signature:
    #     Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu, vs, v0, rs,
    #                      rsarc, s, rminarc, icut, irmax, iradmax, ixmax)
    # f2py sets widtharc as output and inumax optional; call signature is:
    #     fullefflengthsub(kevs, b0, eta2, mu, vs, v0, rs, rsarc, s, rminarc
    #                      icut, irmax, iradmax, ixmax, [inumax])

    fwhms = fm.fullefflengthsub(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s,
                                rminarc, icut, irmax, iradmax, ixmax)

    reserr = fwhms < np.finfo(float).eps * 2
    boxerr = fwhms > (rminarc - np.finfo(float).eps * 2)
    if any(reserr):
        print 'Resolution error! at {}'.format(kevs[np.where(reserr)[0]])
    if any(boxerr):
        print 'Box error! at {}'.format(kevs[np.where(boxerr)[0]])

    return fwhms


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

    v0 = snr.v0()
    rs = snr.rs()
    rsarc = snr.rsarc

    b = snrcat.SYNCH_B
    cm = snrcat.SYNCH_CM
    Cd = snrcat.CD
    nu2 = snrcat.NUKEV * 2.0
    beta = snrcat.BETA
    
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


# ======================================================
# Useful tees -- for code-integrated log file generation
# ======================================================
# (in principle we could even parse the logs to recover model FWHMs if
# something in the code goes wrong -- but we don't want to do that...)

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
