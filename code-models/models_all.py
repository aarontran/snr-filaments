"""
General code to execute model fits for measured FWHMs

(calls f2py wrapper to FullEfflength_mod.f, uses code from Widthfun.py
and inherits code design from earlier prototype `fwhm-process.ipynb`)

Functionality:
    1. manual fitting - let the user twiddle stuff by hand, just a wrapper here,
    no different than calling a.out all the time, just easier
    2. blind fitting - run lmfit with arbitrary initial guesses
    3. grid fitting - generate a grid of values...

Objective:
    1. reproduce SN 1006 numbers with manual fitting
    2. get new SN 1006 numbers with blind fitting
    3. compare blind fits to grid fits (after we're sure code works)
    
    Discuss

    4. Run code on Tycho (need function to parse region files and get shock
    speeds for each location)

Aaron Tran
2014 July 21 (last modified: 2014 July 22)
"""

from __future__ import division

import lmfit
import matplotlib.pyplot as plt
import numpy as np
from numpy import f2py
import scipy as sp
from scipy import optimize
import sys

from fplot import fplot
import snr_catalog as snrcat

# Compile and import Fortran code -- not sure where to do this.
# Could throw this into another short script? remember to call
# fm.readfglists() before anything else though.

#with open('FullEfflength_mod.f', 'r') as f:
#    fsource = ''.join(list(f))
#np.f2py.compile(fsource, modulename='fullmodel', verbose=1)
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


def main():
    """Just testing, for now"""

    kevs = np.array([0.7, 1., 2.])

    flmt = 5
    data, eps = SN1006_DATA[flmt]
    kevs, data, eps = kevs[1:], data[1:], eps[1:]  # Try reproducing SN 1006 numbers

    manual_fit(snrcat.make_SN1006(), kevs, data, eps)

def holding_test():

    # SN 1006 data (how to read this in generally? temporary for testing)
    kevs = np.array([0.7, 1., 2.])
    #data4 = np.array([29, 23.9, 16.6])
    #eps4 = np.array([.9, .39, .45])
    data5 = np.array([33.75, 27.2, 24.75 ])
    eps_data5 = np.array([2.37,.62,.61])

    data, eps = data5, eps_data5
    kevs, data, eps = kevs[1:], data[1:], eps[1:]  # Try reproducing SN 1006 numbers

    # Set SNR physical parameters here
    snr = snrcat.make_SN1006()

    # Search previously generated table for best initial guesses
    # OR, supply your own initial guesses and let er rip
    p = lmfit.Parameters()
    p.add('B0', value=2e-4, min=1e-6, max=1e-2)
    #p.add('eta2', value=1, vary=False)
    p.add('eta2', value=1, min=0, max=1e4)
    p.add('mu', value=1, vary=False)

    # Fit complex model
    res = lmfit.minimize(objectify(width_cont), p,
                   args=(data, eps, kevs, snr),
                   kws={'icut':1},  # Numerical code settings
                   method='leastsq',
                   epsfcn=1e-5)  # this epsfcn is IMPORTANT!!!
                   # Step size set in an ad hoc fashion -- sort of guess and
                   # check.  Goal is to let fit explore parameter space, while 
                   # preventing blowup or immobilization...
        # But, larger steps means longer time to converge, too
        # And, still not guaranteed to find best fit.

    # Print fit information
    lmfit.printfuncs.report_fit(res.params)
    print 'chisqrt = {}'.format(res.chisqr)
    #ci = lmfit.conf_interval(res)  # Only if more than 2 variables?
    #lmfit.printfuncs.report_ci(ci)

    # Fit simple model
    #print 'Simple model'
    #res2 = lmfit.minimize(objectify(width_dump), p,
    #                args=(data, eps, kevs, snr),
    #                method='leastsq')

    # Print fit information
    #lmfit.printfuncs.report_fit(res2.params)
    #ci = lmfit.conf_interval(res2)  # Only if more than 2 variables?
    #lmfit.printfuncs.report_ci(ci)

    # Plot fit with data
    #plot_data()
    en_plt = np.linspace(0.6, 2.1, 10)
    plt.errorbar(kevs, data, eps, fmt='bo')
    plt.plot(en_plt, width_cont(p, en_plt, snr), '-k.')
    plt.xlabel('Observation energy (keV)')
    plt.ylabel('FWHM (arcsec)')
    plt.draw()
    plt.show(block=False)

    # Maybe good to spit out a configuration file here, to log fit inputs and
    # outputs.

    return res


# =============================================================================
# Functions to execute various types of fits
# 1. manual_fit(...) - interactively guess/check parameter values
# 2. tab_eta2(...) and helpers - tabulate parameter space....
# =============================================================================

def fit_from_tab():
    """Fit data from precomputed table of values..."""
    pass

# Here I slaughtered the functions making a mess
# generating different plots and numbers, just exploring parameter space...

def make_tabs():
    """Making tables and plots and things, using grid function"""
    mu_vals = [0, 1/3, 1/2, 2/3, 1, 1.5, 2, 3, 4, 5]
    fmt_vals = [':b', ':g', ':r', ':k', '-k',
                '-r', '-g', '-b', '-c', '-y', '-m']
    #mu_vals = [1]
    #fmt_vals = ['-ok']

    snr = snrcat.make_SN1006()

    for flmt in [1,2,3,4,5]:

        plt.figure(figsize=(7,4))

        kevs = SN1006_KEVS
        data, eps = SN1006_DATA[flmt]
        for mu, fmt in zip(mu_vals, fmt_vals):
            eta2_vals, b0_vals, chisqr_vals = make_tab(snr, kevs, data, eps, mu, fmt)
            plt.plot(eta2_vals, chisqr_vals, fmt, lw=2)

        fplot(r'$\eta_2$ (-)', r'Best fit $\chi^2$', scales=('log','log'))
        plt.title((r'Filament {:d}, simple fits at various $\eta_2$ '
                   '($B_0$ free)').format(flmt))
        plt.legend( tuple(r'$\mu = {:0.2f}$'.format(mu) for mu in mu_vals),
                    loc='best')

        #plt.savefig('sn1006-flmt_{:d}-simple_fit-eta2_space.pdf'.format(flmt))
        plt.draw()
        plt.show()
        plt.clf()


def tab_eta2(snr, kevs, data, eps, mu, fmt):
    """SNR dependent tabulation of FWHM values for fitting / to explore
    parameter space, given some energy bands
    Carve out a swath of points in the eta2-B0 plane for fixed mu.

    First, grid over eta2.
    For each eta2, get range of B0 values that give reasonable FWHMs.
    With B0 values, also save the model FWHMs.

    Verify that results are independent of input FWHMs (mins, maxs, etc)
        Generate two tables w/ slightly different FWHM inputs/ranges
        Fit procedure should give results consistent within error

    Think about full model fit errors?  Once we have best fit B0, eta2, mu for a
    given region/filament, then deal with error.  Here we're just making
    a table of possible/candidate values.
    """

    # Can go ahead and grid over mu in here as well...

    eta2_vals = np.logspace(-2, 3, 20, base=10)  # Must be sorted
    eta2_dict = {}

    p = lmfit.Parameters()
    p.add('mu', value=mu, vary=False)

    for eta2 in eta2_vals:
        p.add('eta2', value=eta2, vary=False)
        p.add('B0', value=150e-6, min = 1e-6, max=1e-2)  # Just initial guess

        # Compute simple fits to "mean" FWHMs, (min+max)/2, to obtain initial B0
        # No weights in fit.
        dmid = (datamax + datamin)/2  # FWHM values in middle of range
        res = lmfit.minimize(objectify(width_dump), pars,
                             args=(dmid, np.ones(len(dmid)), kevs, snr),
                             method='leastsq')

        # List of B0 values, and list of model FWHMs (list of lists)
        eta2_dict[eta2] = tab_eta2_helper_B0(snr, p, kevs, datamin, datamax)

    return eta2_dict


def tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max):
    """
    Supply list of min FWHMs, max FWHMs.  ALL model fit FWHMS must have at
    least one FWHM in the ranges.  If all FWHMs are outside range, then there
    exists a better fit somewhere in parameter space for sure.

    Input:
        pars (lmfit.Parameters): mu, eta2, initial B0 guess, any other parameters needed...

        (alternate approach -- try to grid consistently in B0 and use B0 values
        from previous fit???)

        fwhms_min, fwhms_max must be numpy arrays
    Outputs:
        list of B0 values, list of FWHM value lists (for each B0)
    """

    # TODO: If fits give bad box / bad resolution, it may be worth modifying
    #       FORTRAN code to give 0's or INFs, to ensure they're out of range.
    # If code works (we have good resolution (rminarc small enough),
    # good min/max FWHMs, then we should be okay.

    # TODO: maybe worth logging script output, to check for
    # box/resolution errors.  Shouldn't be a problem though...

    f_minarc = 1.2  # How much margin of safety for rminarc?
    n_min = 20  # Ad hoc guess -- depends on spread of your SNR FWHMs

    # Work in rescaled width space (r in [0,1]), averaged over FWHMs
    def rscale(x):
        return np.mean((x-fwhms_min)/(fwhms_max - fwhms_min))
    def inband(x):
        return all(x < fwhms_max) and all(x > fwhms_min)

    # Get FWHMs from initial guess
    fwhms_init = width_cont(pars, kevs, snr, verbose=False,
                       rminarc = max(fwhms_max)*f_minarc)
    B0_init = pars['B0'].value

    # Reinit B0 if initial guess is bad.  Cannot be TOO bad or this will take a
    # while to run
    # TODO: risk of oscillating forever, back and forth between bad values.
    # (e.g., if B0/1.05 is too small, but B0*1.05 is too big.
    if all(fwhms_init > fwhms_max):
        pars.add('B0', value=B0_init * 1.05)  # min/max not needed, not fitting
        return tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max)
    elif all(fwhms_init < fwhms_min):
        pars.add('B0', value=B0_init / 1.05)
        return tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max)


    # fwhms are in acceptable range.  Start generating B0 values
    r_init = rscale(fwhms_init)

    n_thin = int(np.around(r_init * n_tot))
    n_wide = int(np.around((1 - r_init) * n_tot))
    dr = 1.0 / n_tot

    B0_thin, fwhms_thin = find_pts(r_init, 0, -1*dr)  # Decrease r, increase B0
    B0_wide, fwhms_wide = find_pts(r_init, 1, dr)  # Increase r, decrease B0

    # These must all be regular Python lists
    B0_vals = B0_thin + [B0_init] + B0_wide
    fwhms_all = fwhms_thin + [fwhms_init] + fwhms_wide

    return B0_vals, fwhms_all

"""
To illustrate the general utility of the following functions,
fill_pts(...), mind_the_gap(...), try this code:
(was used for debugging -- making sure it works when iterating backwards)

import numpy as np
import matplotlib.pyplot as plt
import models_all as ma

# Cannot start at (0,0) due to zero derivative
x, y = ma.fill_pts(0.1, 0.1**2, 10, 1, lambda x: x**2, dx_max=1)
xnew, ynew = ma.mind_the_gaps(x, y, 0.5, lambda x: x**2)
plt.plot(x,y,'-ro',ms=10)
plt.plot(xnew,ynew,'*b',ms=8)
plt.show()

x, y = ma.fill_pts(-0.1, (-0.1)**3, -15, 1, lambda x: x**3, dx_max=1)
xnew, ynew = ma.mind_the_gaps(x, y, 0.5, lambda x: x**3)
plt.plot(x,y,'-ro',ms=10)
plt.plot(xnew,ynew,'*b',ms=8)
xm = np.linspace(1, -(15)**(1./3), 100)
plt.plot(xm, xm**3, '-g')
plt.show()
"""


def fill_pts(x0, y0, y_final, dy_step, f, dx_max=float('inf'),
             epsfcn_init=1e-2):
    """Find pts (x,y) with y-values spanning interval [y0, y_final]

    Just Newton-Raphson iteration, but keeping intermediate values.
    dy_step is just a guideline for stepping torwards y_final
    to enforce data spacing, use mind_the_gaps

    WARNINGS: function must be monotonic, well-behaved, and have no extrema
    (i.e., derivative strictly nonzero and single-signed).
    There is no upper bound on the number of points generated.
    Code developed specifically for Sean's model fitting.

    Inputs:
    Outputs:
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
    x_prev, y_prev = x_curr, y_curr = x0, y0  # Just for starters

    while sgn * y_curr < sgn * y_final:  # ad hoc stopping threshold
        # if dy_step > 0, proceed if y_curr < y_final
        # if dy_step < 0, proceed if y_curr > y_final
        if len(x_vals) > 1:
            jcbn = (y_curr - y_prev) / (x_curr - x_prev)  # Approx first deriv
        else:  # Initialize first deriv. if only (x0, y0) known
            dx_init = epsfcn_init * abs(x0)
            jcbn = (f(x0 + dx_init) - y0)/dx_init  # Could store width measurement...

        dx = next_step(jcbn)

        # Shift last iteration values back
        x_prev, y_prev = x_curr, y_curr
        # Add new data point
        x_curr = x_prev + dx
        y_curr = f(x_curr)  # Store this width measurement
        x_vals = np.append(x_vals, x_curr)
        y_vals = np.append(y_vals, y_curr)

    return x_vals, y_vals


def mind_the_gaps(x, y, dy_max, f):
    """Fill data until y-spacing is always less than dy_max
    Assumptions: x, y are sorted; y = f(x)

    Whenever a gap (dy > dy_max) is found, iteratively sample f at x-coordinate
    midpoints until all gaps are smaller than dy_max

    Input:
        x, y (np.array): data values to fill in
        dy_max (float): maximum y spacing, must be greater than zero
        f (function): function with call signature f(x).  Expect y = f(x)
    Output:
        x, y (np.array) with new data values
    """
    dy = np.abs(np.diff(y))
    while any(dy > dy_max):
        mask = dy > dy_max 
        indgaps = np.where(mask)[0]

        x_midpts = (x[1:][mask] + x[:-1][mask]) / 2
        y_midpts = np.array([f(xm) for xm in x_midpts])  # Store widths here

        x = np.insert(x, indgaps + 1, x_midpts)
        y = np.insert(y, indgaps + 1, y_midpts)
        dy = np.abs(np.diff(y))

    return x, y


def simple_fits(snr, kevs, data, eps, mu):
    """Fit both eta2 and B0 together, simple model only.
    Should reproduce Table 7 of Ressler!
    """

    p.add('B0', value=150e-6, min=1e-6, max=1e-2)
    p.add('eta2', value=1, min=0.01, max=1000)
    res = lmfit.minimize(objectify(width_dump), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq')
    print 'Best fit values with both eta2, B0 free'
    print 'eta2 = {:.2f}, B0 = {:.2f}, mu = {:.2f}'.format(p['eta2'].value,
            p['B0'].value * 1e6, p['mu'].value)
    print 'measured: {}'.format(data)
    print ('simple fit: {}; '.format(width_dump(p, kevs, snr))
            + 'chisqr = {:.2f}'.format(res.chisqr))

    return res


def manual_fit(snr, kevs, w_meas, w_eps):
    """Interactively prompts user for values of B0, eta2, mu
    Prints FWHMs, residuals, and chi-squared values in the process
    Also prints m_E, computed point to point (just for comparison)

    Note: doesn't really deal with bad input

    Input: snr object to be fitted
    Output: None? or, best attempted fit parameters to date
    """

    print 'Manual fitting routine: enter q to print best fit and quit.\n'

    chisq_best = float('inf')
    w_model_best = None
    p_best = None

    while True:
        B0 = get_float('Enter B0 (G): ')
        eta2 = get_float('Enter eta2 (-): ')
        mu = get_float('Enter mu (-): ')

        if any(np.isnan([B0, eta2, mu])):
            break  # Received exit/quit request

        p = lmfit.Parameters()
        p.add('B0', value=B0)
        p.add('eta2', value=eta2)
        p.add('mu', value=mu)

        adj = get_float('Change model settings? Enter 0/1 for no/yes: ')
        if adj == 0:
            w_model = width_cont(p, kevs, snr)
        else:
            rminarc = get_float('rminarc (default: 60): ')
            icut = get_float('icut (default: 1): ')
            #irmax = get_float('irmax (default: ): ')
            #iradmax = get_float('iradmax (default: 60): ')
            #ixmax = get_float('ixmax (default: 60): ')
            w_model = width_cont(p, kevs, snr, rminarc=rminarc, icut=icut)

        chisq = chi_squared(w_meas, w_eps, w_model)
        wresid = (w_model - w_meas)/w_eps

        print ('\nObsvd fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*w_meas)
        print ('Model fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*w_model)
        print ('Wghtd resid: ' + len(kevs)*'{:0.3f}, ')[:-2].format(*wresid)
        print 'Chi^2: {:0.3f}'.format(chisq)

        plt.clf()
        plt.errorbar(kevs, w_meas, w_eps, fmt='bo')
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
            print ('Best fwhms:  ' + len(kevs)*'{:0.2f}, ')[:-2].format(*w_model_best)
            print 'Chi^2: {:0.3f}'.format(chisq_best)

        print '--------------------------------'
    
    if p_best is not None:
        print '\nDone with manual fit. Best fit parameters are:'
        for key in p_best.keys():
            print '{} = {:g}'.format(p_best[key].name, p_best[key].value)


def chi_squared(y, y_err, y_model):
    """Compute chi-squared statistic"""
    return np.sum( ((y_model - y) / y_err)**2 )


def get_float(prompt):
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
    v0 = snr.v0
    rs = snr.rs
    rsarc = snr.rsarc
    s = snr.s

    if rminarc is None:
        rminarc = snr.rminarc

    if verbose:
        print ('\tFunction call with B0 = {:0.3f} muG; eta2 = {:0.3f}; '
               'mu = {:0.3f}').format(B0*1e6, eta2, mu)

    # FORTRAN subroutine call signature:
    #     Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu, vs, v0, rs,
    #                      rsarc, s, rminarc, icut, irmax, iradmax, ixmax)
    # f2py sets widtharc as output and inumax optional; call signature is:
    #     fullefflengthsub(kevs, b0, eta2, mu, vs, v0, rs, rsarc, s, rminarc
    #                      icut, irmax, iradmax, ixmax, [inumax])

    return fm.fullefflengthsub(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s,
                               rminarc, icut, irmax, iradmax, ixmax)


def width_dump(params, kevs, snr):
    """Width function from catastrophic dump model (equation 6), code taken
    from Sean's code Widthfun.py
    
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

    v0 = snr.v0
    rs = snr.rs
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
    width = beta * 2*D/v0 / (np.sqrt(1 + (4*D/v0**2/tsynch)) - 1)
    
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


if __name__ == '__main__':
    main()
