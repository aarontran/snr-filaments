"""
utilities for
1. finding FWHMs,
2. computing FWHM errors by stretching...

For each profile (every band, every region), we need to:
1. calculate FWHM
2. estimate uncertainties on FWHM
3. check FWHM errors are reasonable (i.e., $\chi^2$ space is not too weird)
4. give information to allow plotting/visualization of uncertainties.

I am iffy about the quality of some of these functions -- some can probably be
cleaned / simplified, for clarity

Aaron Tran
July 2014 (moved to .py, Sept 2014)
"""

import numpy as np
import scipy as sp
from scipy import optimize

from crvfit import chi2, chi2red

# mega method to get everything at once because I'm lazy

def get_fwhm_all(x, y, y_err, f, pars, want_err=True,y_cap=None,dx=5,**kwargs):
    """After running a fit routine, compute useful numbers associated with FWHM

    Input:
        x, y, y_err: measured profile data
        f, pars: profile fit function, best fit parameters
        want_err (boolean): compute errors or not?
        y_cap (float): artificially change max for FWHM calculation
        dx (float): search range (arcsec) for function max in nbhd of data max
                    (5 arcsec appears to work well for Tycho)
        Other **kwargs passed to stretch computation functions
    Output:
        fwhm, (fwhm_a, fwhm_b), (neg_err, pos_err)
    """
    x0 = x[np.argmax(y)]
    xmax = find_fmax(x0-dx, x0+dx, f, pars)  # Search within +/- dx arcsec. of x0
    # This is the FUNCTION's xmax, not the DATA xmax

    fwhm, fwhm_a, fwhm_b = get_fwhm(np.amin(x), np.amax(x), f, pars,
                                    x_max=xmax, y_cap=y_cap)

    if want_err and np.isfinite(fwhm):
        xi_neg, xi_pos = get_stretch_limits(x, y, y_err, f, pars, xmax, **kwargs)
        fwhm_errs = get_fwhm_errs(fwhm, fwhm_a, fwhm_b, xi_neg, xi_pos, xmax,
                                  **kwargs)
    else:
        fwhm_errs = (float('nan'), float('nan'))

    return fwhm, (fwhm_a, fwhm_b), fwhm_errs

# FWHM calculation functions

def find_fmax(a, b, f, pars):
    """LOCAL maximum of a scalar function, within a given interval"""
    res = sp.optimize.minimize_scalar(lambda x: -1*f(x, *pars), bounds=(a, b),
                                      method='Bounded')
    return res.x

def get_fwhm(a, b, f, pars, x_max=None, y_cap=None):
    """Calculate FWHM given bracketing interval, maximum, and function.

    Assumptions -- there exists a single local maximum, from which the function
    MONOTONICALLY DECREASES on both sides

    If y_max is supplied, use THAT instead of actual maximum!
    Then it's no longer a true FWHM for the function, but may be better
    estimate of FWHM for the data

    Inputs:
        a, b (float): bounding interval
        f (function): scalar function in first argument, with form f(x, *popt)
        pars (list): list of parameters for function f
        x_max (float): location of maximum if known (skip maximization step)
        y_cap (float): fixed "maximum" value, to get full-width-at-half-"y_cap" instead

    Output:
        3-tuple of FWHM and left, right locations of FWHM measurement
    """
    if x_max is None:
        x_max = find_fmax(a, b, f, pars)

    if y_cap is None:
        y_hm = f(x_max, *pars) / 2.
    else:
        y_hm = y_cap / 2.

    y_lims = f(np.array([a,b]), *pars)  # Assumes that f is vectorized in x
    if any(y_lims >= y_hm):  # Function must fall below half max on both sides
        print 'WARNING: cannot calculate FWHM, peak not resolved in interval'
        return float('NaN'), float('NaN'), float('NaN')

    fwhm_left = sp.optimize.brentq(lambda x: f(x, *pars) - y_hm, a, x_max)
    fwhm_right = sp.optimize.brentq(lambda x: f(x, *pars) - y_hm, x_max, b)
    return fwhm_right - fwhm_left, fwhm_left, fwhm_right

def get_fwhm_errs(fwhm, fwhm_left, fwhm_right, xi_neg, xi_pos, x0, **kwargs):
    """FWHM limits based on stretch limits (in x-coordinate only)
    Applies an inverse stretch operation to transform
    old coordinates to new coordinates

    E.g., f_stretch(new_fwhm_b) = f(x'(new_fwhm_b)) = f(old_fwhm_b)
    So we solve the inverse problem old_fwhm_b = x'(new_fwhm_b)

    Inputs:
        fwhm: fwhm
        fwhm_a, fwhm_b: x coordinate bounds on fwhm location
        xi_neg, xi_pos: some confidence limits on xi
        x0: center for xi-stretch
        **kwargs: sent to inv_stretch_x (set xconst)
    Output:
        errors on FWHM
    """
    intrv = np.array([fwhm_left, fwhm_right])

    min_intrv = inv_stretch_x(intrv, xi_pos, x0, **kwargs)
    fwhm_min = min_intrv[1] - min_intrv[0]

    max_intrv = inv_stretch_x(intrv, xi_neg, x0, **kwargs)
    fwhm_max = max_intrv[1] - max_intrv[0]

    err_neg = fwhm - fwhm_min
    err_pos = fwhm_max - fwhm

    return err_neg, err_pos

# Main work here -- stretch function until confidence limit reached

def get_stretch_limits(x, y, y_err, f, pars, xmax, **kwargs):
    """Limits on stretch parameter xi for 90% confidence interval
    x-coordinate stretch function follows Ressler et al.
    pars specifies least-squares best fit parameters
    **kwargs sent to safe_stretch, stretch_func
    """

    def chi2_threshold(xi):
        """Threshold for chi^2 uncertainty. Output > 0 if \Delta\chi^2 > 2.7"""
        def chi2_stretch(xi):
            return chi2(x, y, y_err, stretch_func(f, xi, xmax, **kwargs), pars)
        chi2_stretch = np.vectorize(chi2_stretch)

        chi2_min = chi2(x, y, y_err, f, pars)  # same as chi2_stretch(xi=0)

        return chi2_stretch(xi) - chi2_min - 2.7  # Magic number...

    # Increasing positive xi (decreasing FWHM) uncertainty
    # Empirically, positive xi is better behaved than negative xi
    xi_max = 0.01
    while chi2_threshold(xi_max) < 0:
        xi_max += 0.01
    xi_bound_pos = sp.optimize.brentq(chi2_threshold, 0, xi_max)
    if not safe_stretch(xi_bound_pos, xmax, **kwargs):
        print 'WARNING: xi may be too large, pathological behavior may occur for small x'

    # Decreasing negative xi (increasing FWHM) uncertainty
    # Empirically, this can blow up quickly (stretch becomes meaningless if too far)
    xi_min = -0.01
    while chi2_threshold(xi_min) < 0:
        xi_min -= 0.01
    xi_bound_neg = sp.optimize.brentq(chi2_threshold, xi_min, 0)
    if not safe_stretch(xi_bound_neg, xmax, **kwargs):
        raise Exception('ERROR: xi negative and too large, pathological behavior WILL OCCUR (flipping)')

    return xi_bound_neg, xi_bound_pos

# Stretching functions to compute errors

def stretch_func(f, xi, x_max, **kwargs):
    """Stretched version of a function with call signature f(x, *pars)
    If function f is vectorized in x, output function will also be vectorized"""
    return lambda x, *pars: f(stretch_x(x, xi, x_max, **kwargs), *pars)

def safe_stretch(xi, x0, xconst=50):
    """Check if xi is in the safe range
    I.e., the range in which
    1. stretch works as desired on interval [0, xconst]
    2. stretch is 1-to-1 on said interval
    """
    if x0 > xconst:
        raise Exception('Invalid stretch x0={} > xconst'.format(x0))
    return abs(xi) < (xconst - x0)/x0

def stretch_x(x, xi, x0, xconst=50):
    """Parabolic stretch centered on x=x0, fixed points at x=0, x=x0
    Follows Ressler et al. [2014]

    Constant here is 50 arcsec instead of 200 arcsec.
    Doesn't matter, mainly scales the input range of xi

    No error checking -- you put bad xi, you get bad stretch

    Inputs:
        x (float / iterable): original x-values
        xi (float): stretch factor
        x0 (float): fixed point of stretch
    Output:
        stretched version of x (float/iterable depending on input)
    """
    return x * (1 + xi * (x - x0)/(xconst - x0))

def inv_stretch_x(xs, xi, x0, xconst=50):
    """Calculates the desired value of x, given a stretched x
    Returns the "sensible" value for our stretch

    If a=xi >0, then positive root is bigger, positive (2nd root)
    If a=xi <0, then positive root is smaller (1st root)

    Can return a NaN value if stretch cannot be inverted
    (physical meaning: stretch breaks down at large distances)
    Vectorized in xs, by nature (but not in other arguments)
    """
    a = xi
    b = xconst - x0*(1 + xi)
    c = -(xconst - x0) * xs
    d = b**2 - 4*a*c
    posroot = (-b + np.sqrt(d)) / (2*a)

    if not safe_stretch(xi, x0, xconst):
        if xi < 0:
            raise Exception('xi={} < {}'.format(xi))
        else:
            print 'Warning: possibly bad inverse stretch (xi too large)'

    if xi == 0:
        return xs
    return posroot

if __name__ == '__main__':
    pass
