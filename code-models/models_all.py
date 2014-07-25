"""
General code to execute model fits for measured FWHMs

(calls f2py wrapper to FullEfflength_mod.f, uses code from Widthfun.py
and inherits code design from earlier prototype `fwhm-process.ipynb`)

Functionality:
    1. manual fitting - let the user twiddle stuff by hand, just a wrapper.
       no different than calling a.out all the time, just easier
    2. blind fitting - run lmfit with arbitrary initial guesses.
       this is a BAD IDEA on anything other than simple model.
    3. grid fitting - generate a grid of values...

Objective:
    1. reproduce SN 1006 numbers with manual fitting
    2. get new SN 1006 numbers with blind fitting
    3. compare blind fits to grid fits (after we're sure code works)

    Discuss... make tables.

    4. Run code on Tycho (need function to parse region files and get shock
       speeds for each location... can throw that function in here or something.
       Or, separate execution from implementation -- iterate over
       SNRs, filaments, regions in a separate file?

Note: the fitting formalism used here (tabulating function values over a
subset of parameter space by careful gridding + good initial guess), might be
useful in other problems?

Aaron Tran
2014 July 21 (last modified: 2014 July 25)
"""

from __future__ import division

import cPickle as pickle
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

    flmt = 5
    kevs = SN1006_KEVS
    data, eps = SN1006_DATA[flmt]
    kevs, data, eps = kevs[1:], data[1:], eps[1:]  # Try reproducing SN 1006 numbers

    manual_fit(snrcat.make_SN1006(), kevs, data, eps)


def holding_test():
    """Testing some fits -- just messing around with numbers"""
    flmt = 5
    kevs = SN1006_KEVS
    data, eps = SN1006_DATA[flmt]

    kevs, data, eps = kevs[1:], data[1:], eps[1:]  # Try reproducing SN 1006 numbers

    # Set SNR physical parameters here
    snr = snrcat.make_SN1006()

    # Search previously generated table for best initial guesses
    # OR, supply initial guesses
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


def simple_fits(snr, kevs, data, eps, mu):
    """Fit both eta2 and B0 together, simple model only.
    Should reproduce Table 7 of Ressler!
    """
    # this is strictly temporary -- to be adapted into a larger routine

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


# ===============================
# Different fitting procedures
# manual_fit(...), table_fit(...)
# ===============================

def table_fit():
    """Fit data from precomputed table of values..."""
    pass



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


# =========================================
# Tabulate FWHM values from full model code
# =========================================

def make_tabs():
    """Making tables and plots and things, using grid function.
    Temporary -- previously used for some tests/plots.
    Will iterate over mu and call gridding functions when done
    """
    pass # superseded by tab_eta2()? 


def test_gen_table():
    """Try generating a table for SN 1006, and save
    ONLY RUN THIS ONCE LEST YOU OVERWRITE..."""
    tab = tab_eta2()
    

    try:
        fname = 'sn1006_grid-3-40-20_2014-07-25.pkl'
        with open(fname, 'w') as fpkl:
            pickle.dump(tab, fpkl)
    except:
        print 'dammit, pickle dump failed'

    vistable(tab)
    return tab  # TO BE PARANOID


def vistable(tab):
    """Make plot of grid, using tab output from tab_eta2()"""
    mu_vals = tab.keys()
    mu_dict = tab[mu_vals[0]]
    for eta2 in mu_dict.keys():
        B0_vals = mu_dict[eta2][0]  # [1] for FWHM values!
        plt.plot(eta2*np.ones(len(B0_vals)), B0_vals, 'ob')

    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$\eta_2$'); plt.ylabel(r'$B_0$')
    plt.show()


def tab_eta2():
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

    #mu_vals = [0, 1/3, 1/2, 2/3, 1, 1.5, 2, 3, 4, 5]
    #eta2_vals = np.logspace(-2, 3, 20, base=10)  # Must be sorted

    #mu_vals = [0, 1/3, 1/2, 1, 1.5, 2]  # Following Sean
    mu_vals = [1/3, 1/2, 1]  # Kolmogorov, Kraichnan, Bohm
    eta2_vals = np.logspace(-2, 2, 40, base=10)
    n_B0 = 20

    # Estimate: this 3*40*20 * 3sec/call / 3600 = ~ 2 hours to run

    snr = snrcat.make_SN1006()
    # Generate fwhm_range
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    #eps = np.array([SN1006_DATA[flmt][1] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data, axis=0)
    data_max = np.amax(data, axis=0)
    #eps = np.mean(eps, axis=0)
    
    mu_dict = {}

    print 'Gridding with FWHM limits:'
    print 'min: {}'.format(data_min)
    print 'max: {}'.format(data_max)

    for mu in mu_vals:
        print 'Starting grid for mu = {}'.format(mu)

        p = lmfit.Parameters()
        p.add('mu', value=mu, vary=False)

        eta2_dict = {}
        for eta2 in eta2_vals:
            print 'Adding points for eta2 = {}'.format(eta2)

            p.add('eta2', value=eta2, vary=False)
            p.add('B0', value=150e-6, min = 1e-6, max=1e-2)  # Just initial guess

            # Compute simple fits to "mean" FWHMs, (min+max)/2, to obtain initial B0
            # No weights in fit.
            dmid = (data_max + data_min)/2  # FWHM values in middle of range
            res = lmfit.minimize(objectify(width_dump), p,
                                 args=(dmid, np.ones(len(dmid)), kevs, snr),
                                 method='leastsq')

            # List of B0 values, and list of model FWHMs (list of lists) (np.array)
            eta2_dict[eta2] = tab_eta2_helper_B0(snr, p, kevs, data_min,
                                                 data_max, n_B0)

        mu_dict[mu] = eta2_dict

    # write final table to pickle, just temporarily

    return mu_dict


def tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot):
    """
    Supply list of min FWHMs, max FWHMs (rough guidelines for gridding).  ALL model fit FWHMS must have at
    least one FWHM in the ranges (note no longer satisfied...)
    If all FWHMs are outside range, then there
    exists a better fit somewhere in parameter space for sure.

    Procedure:
    1. Improve initial guess for B0 using simple model fit to "average" fwhms
       (midway between fwhms_min, fwhms_max).
    2. From new initial guess for B0, generate grid over many values of B0.
       Store values of B0 and FWHMs from complex model.

    TODO: If fits give bad box / bad resolution, it may be worth modifying
          FORTRAN code to give 0's or INFs, to ensure they're out of range.
    If code works (we have good resolution (rminarc small enough),
    good min/max FWHMs, then we should be okay.

    TODO: maybe worth logging script output, to check for
    box/resolution errors.  Shouldn't be a problem though...

    TODO: rminarc is fixed at f_minarc*max(fwhms_max)...
    should be good enough for our work, BUT if there is a case of wide spread
    in SNR FWHMs or otherwise this may be a problem.

    Input:
        pars (lmfit.Parameters): mu, eta2, initial B0 guess, any other parameters needed...

        (alternate approach -- try to grid consistently in B0 and use B0 values
        from previous fit???)

        fwhms_min, fwhms_max must be numpy arrays
    Outputs:
        list of B0 values, list of FWHM value lists (for each B0)
    """

    f_minarc = 1.2  # How much margin of safety for rminarc?
    n_min = 20  # Ad hoc guess -- depends on spread of your SNR FWHMs

    # Work in rescaled width space (r in [0,1]), averaged over FWHMs
    def rscale(x):
        return np.mean((x-fwhms_min)/(fwhms_max - fwhms_min))
    def f_rscale(grid_B0):
        """Simplify full model function for grid fitting"""
        pars.add('B0', value=grid_B0)  # Vary B0, other parameters same
        fwhms = width_cont(pars, kevs, snr,
                           rminarc = max(fwhms_max)*f_minarc)
        #fwhms = width_dump(pars, kevs, snr)  # FOR FAST TESTING!
        return rscale(fwhms), fwhms
    #def inband(x):
    #    return all(x < fwhms_max) and all(x > fwhms_min)

    # Get FWHMs from initial guess, reinitialize B0 if initial guess is bad
    r_init, fwhms_init = f_rscale(pars['B0'].value)  # Reinit B0 guess
    B0_init = pars['B0'].value
    if all(fwhms_init > fwhms_max):  # WARNING: risk of oscillating forever
        pars.add('B0', value=B0_init * 1.05)  # min/max not needed, not fitting
        return tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max)
    elif all(fwhms_init < fwhms_min):
        pars.add('B0', value=B0_init / 1.05)
        return tab_eta2_helper_B0(snr, pars, kevs, fwhms_min, fwhms_max)

    print 'Using initial B0 value {}'.format(B0_init)
    print 'Now finding B0 values with FWHMs in range.'

    # fwhms are in acceptable range.  Start generating B0 values
    # Idea is to grid evenly over rescaling of FWHMS (rscale)

    n_thin = int(np.around(r_init * n_tot))
    n_wide = int(np.around((1 - r_init) * n_tot))
    dr = 1.0 / n_tot  # y-coord for fill_pts

    # Calculate values
    B0_wide, r_wide, fwhms_wide = fill_pts(B0_init, r_init, fwhms_init, 1, dr, f_rscale,
                                           dx_max=10e-6)  # Increase r, decrease B0
    B0_thin, r_thin, fwhms_thin = fill_pts(B0_init, r_init, fwhms_init, 0, dr, f_rscale,
                                           dx_max=10e-6)  # Decrease r, increase B0

    print 'Filling in FWHM range to achieve desired spacing'

    # Combine lists, remove double included B0_init/r_init and sort
    B0_all = np.append(B0_thin, B0_wide)
    fwhms_all = np.append(fwhms_thin, fwhms_wide, axis=0)
    r_all = np.append(r_thin, r_wide)


    B0_all, uniqsrt = np.unique(B0_all, return_index=True)
    fwhms_all = fwhms_all[uniqsrt]
    r_all = r_all[uniqsrt]

    # Fill in the gaps
    B0_all, r_all, fwhms_all = mind_the_gaps(B0_all, r_all, fwhms_all, dr, f_rscale)

    return B0_all, fwhms_all  # don't care about r values after we're done


def fill_pts(x0, y0, yaux0, y_final, dy_step, f, dx_max=float('inf'),
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
        may not be sorted -- you deal with this.
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
            jcbn = (f(x0 + dx_init)[0] - y0)/dx_init  # Could store width measurement...

        dx = next_step(jcbn)

        # Shift last iteration values back
        x_prev, y_prev = x_curr, y_curr

        # Add new data point
        x_curr = x_prev + dx
        y_curr, yaux_curr = f(x_curr)  # Store this width measurement

        x_vals = np.append(x_vals, x_curr)
        y_vals = np.append(y_vals, y_curr)
        yaux_vals = np.append(yaux_vals, [yaux_curr], axis=0)

    return x_vals, y_vals, yaux_vals


def mind_the_gaps(x, y, yaux, dy_max, f):
    """Fill data until y-spacing is always less than dy_max
    Assumptions: x, y are sorted; y = f(x)

    Whenever a gap (dy > dy_max) is found, iteratively sample f at x-coordinate
    midpoints until all gaps are smaller than dy_max

    Input:
        x, y (np.array): data values to fill in
        dy_max (float): maximum y spacing, must be greater than zero
        f (function): function with call signature f(x).  Expect y, yaux = f(x)
    Output:
        x, y (np.array) with new data values

        yaux = fwhms
    """
    dy = np.abs(np.diff(y))
    while any(dy > dy_max):
        mask = dy > dy_max 
        indgaps = np.where(mask)[0]

        x_midpts = (x[1:][mask] + x[:-1][mask]) / 2

        # Parse main + auxiliary output
        ytup_out = [f(xm) for xm in x_midpts]
        y_midpts = np.array([tup[0] for tup in ytup_out])
        yaux_midpts = np.array([tup[1] for tup in ytup_out])
        # y_midpts = np.array([f(xm) for xm in x_midpts])  # Store widths here

        x = np.insert(x, indgaps + 1, x_midpts)
        y = np.insert(y, indgaps + 1, y_midpts)
        yaux = np.insert(yaux, indgaps + 1, yaux_midpts, axis=0)
        dy = np.abs(np.diff(y))

    return x, y, yaux


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
