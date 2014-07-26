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


def main():
    """Just for now"""
    generate_SN1006_table()


# ======================================================
# Code to execute fits (this belongs in a separate file)
# ======================================================

def run_table_fit():
    """Load thing"""

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS

    with open('tables/sn1006_grid-3-40-15_2014-07-25.pkl', 'r') as fpkl:
        tab = pickle.load(fpkl)

    for flmt in [1]: #[1,2,3,4,5]:
        print 'Filament {}'.format(flmt)
        print 'mu\teta2\t\t\tB0'
        data, eps = SN1006_DATA[flmt]

        for mu in [1/3, 1/2, 1]:
            res = table_fit(snr, kevs, data, eps, mu, tab)
        plt.show()

            #print '{:0.2f}\t{:0.2f} +/- {:0.2f}\t\t{:0.2f} +/- {:0.2f}'.format(
            #        res.values['mu'],
            #        res.values['eta2'], res.params['eta2'].stderr,
            #        res.values['B0']*1e6, res.params['B0'].stderr*1e6)

            # Next, error from confidence intervals...
            # Error bars on B0 and eta2 will be tremendous in many cases.


def check_effects():
    """Make numbers blagh
    What happens when we change
    compression ratio (i.e. v0) ?
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
            

# ============================================
# Fit from just initial guesses (simple model)
# ============================================

# May be good to somehow save configuration files --... logt settings, SNR
# parameters, et cetera...

def simple_fit(snr, kevs, data, eps, mu):
    """Fit both eta2 and B0 together, simple model only.
    Should reproduce Table 7 of Ressler!
    """
    p = lmfit.Parameters()
    p.add('mu', value=mu, vary=False)
    p.add('B0', value=150e-6, min=1e-6, max=1e-2)
    p.add('eta2', value=1, min=0.01, max=1000)
    res = lmfit.minimize(objectify(width_dump), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq')

    #print 'Best fit values with both eta2, B0 free'
    #print 'eta2 = {:.2f}, B0 = {:.2f}, mu = {:.2f}'.format(p['eta2'].value,
    #        p['B0'].value * 1e6, p['mu'].value)
    #print 'measured: {}'.format(data)
    #print ('simple fit: {}; '.format(width_dump(p, kevs, snr))
    #        + 'chisqr = {:.2f}'.format(res.chisqr))

    return res


# ======================================================
# Fit, optimizing with precached table (full model code)
# ======================================================

def table_fit(snr, kevs, data, eps, mu, tab):
    """Fit data from precomputed table of values...
    
    mu must be an element of tab.keys(); we neglect fit in mu space.
    tab should be opened prior
    """

    # Find best fit values in table, eta2 and B0 for each mu or whatever
    # The fear is that, if we are not sampling finely enough in B0 or eta2
    # We will miss the correct global minimum...
    eta2_dict = tab[mu]

    eta2_bestchisqrs = np.array([])
    eta2_bestB0 = np.array([])

    eta2_vals = np.sort(eta2_dict.keys())
    print eta2_vals

    for eta2 in eta2_vals:
        B0_vals, fwhms = eta2_dict[eta2]  # Unpack tuple

        chisqr_vals = map(lambda x: chi_squared(data, eps, x), fwhms)
        ind = np.argmin(chisqr_vals)

        #eta2_bestchisqrs = np.append(eta2_bestchisqrs, chisqr_vals[ind])
        #eta2_bestB0 = np.append(eta2_bestB0, B0_vals[ind])

        p = lmfit.Parameters()
        p.add('mu', value=mu, vary=False)
        p.add('B0', value=B0_vals[ind], min=1e-6, max=1e-2)
        p.add('eta2', value=eta2, vary=False)
        res = lmfit.minimize(objectify(width_cont), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq',
                             epsfcn=1e-6)  # Try messing with epsfcn...

        eta2_bestchisqrs = np.append(eta2_bestchisqrs, res.chisqr)
        eta2_bestB0 = np.append(eta2_bestB0, p['B0'].value)

        print ('Table best fit: B0 = {:0.2f}, chisqr = '
               '{:0.2f}').format(B0_vals[ind]*1e6, chisqr_vals[ind])
        print ('Fitted best fit: B0 = {:0.2f}, chisqr = '
               '{:0.2f}').format(p['B0'].value*1e6, res.chisqr)

    fmt = {1/3: '-or', 1/2:'-og', 1:'-ob'}
    plt.plot(eta2_vals, eta2_bestchisqrs, fmt[mu])
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel(r'$\eta_2$'); plt.ylabel(r'$\chi^2$')
    #plt.show()

    # From best fit table values, run one more fit

    #res = lmfit.minimize(objectify(width_cont), p,
    #               args=(data, eps, kevs, snr),
    #               kws={'icut':1},  # Numerical code settings
    #               method='leastsq',
    #               epsfcn=1e-5)  # Try messing with epsfcn...

    # Make plots of chi-squared as a function of eta2
    # Give fit information, estimate error (or let lmfit try to do so)
    #lmfit.printfuncs.report_fit(res.params)
    #print 'chisqrt = {}'.format(res.chisqr)
    #ci = lmfit.conf_interval(res)  # Only if more than 2 variables?
    #lmfit.printfuncs.report_ci(ci)
    # return res


# ============================================
# Interactive manual fitting (full model code)
# ============================================

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

def merge_tables(tab1, tab2):
    """does what it says"""
    pass


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


def generate_SN1006_table():
    """Try generating a table for SN 1006, and save
    TODO: Once done -- you should move table to /tables and set write only...
    code to merge table pickles coming sometime...
    TODO: would also be good to save data in plaintext form, somehow
    (4+ columns of mu, eta2, B0, fwhm values)
    TODO: no more feature creep dammit.
    """

    mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
    #mu_vals = np.array([1./3, 1./2, 1])  # Kolmogorov, Kraichnan, Bohm
    eta2_vals = np.logspace(-2, 2, 50, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data/1.51, axis=0)
    data_max = np.amax(data*1.51, axis=0)

    fname = 'sn1006_grid-6-100-20_2014-07-25.pkl'

    tab = maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0, fname)

    #vistable(tab)
    return tab  # TO BE PARANOID


def maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals, n_B0, fname_out):
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

    # Start logging to files (search for errors...)
    # Beware -- this messes with the stdout improperly...
    # doesn't return control correctly if ctrl-c'd.
    log_fname = fname_out + '.log'
    errlog_fname = fname_out + '.errlog'  # temporary.... use os.path
    stdout = TeeStdout(log_fname, 'w')
    stderr = TeeStderr(errlog_fname, 'w')

    np.set_printoptions(precision=2)
    print '\nTabulating full model code FWHMs for {}'.format(snr.name)
    print 'Started: {}'.format(datetime.now())
    print 'Resolution in mu, eta2, B0: {}, {}, {}+'.format(len(mu_vals),
            len(eta2_vals), n_B0)
    print 'Mu values are {}'.format(mu_vals)
    print 'eta2 values are {}'.format(eta2_vals)
    print 'Gridding with FWHM limits:'
    print 'min: {}'.format(data_min)
    print 'max: {}'.format(data_max)

    mu_dict = {}
    for mu in mu_vals:
        p = lmfit.Parameters()
        p.add('mu', value=mu, vary=False)

        mu_dict[mu] = eta2_dict = {}
        for eta2 in eta2_vals:
            print '\nAdding points for mu = {:0.2f}, eta2 = {:0.2f}'.format(mu, eta2)
            print '-------------------------------------------'

            p.add('eta2', value=eta2, vary=False)
            p.add('B0', value=150e-6, min = 1e-6, max=1e-2)

            # Give initial guess for B0, using simple model fit
            dmid = (data_max + data_min)/2
            res = lmfit.minimize(objectify(width_dump), p, method='leastsq',
                                 args=(dmid, np.ones(len(dmid)), kevs, snr))

            # Do the heavy work of computing FWHMs for many B0 values
            eta2_dict[eta2] = maketab_gridB0(snr, p, kevs, data_min, data_max,
                                             n_B0)

            # Save data, but don't interrupt tabulation if fails
            try:
                with open(fname_out, 'w') as fpkl:
                    pickle.dump(mu_dict, fpkl)
            except:
                print 'ERROR: could not save table to file'

    # Stop logger
    print 'Finished: {}'.format(datetime.now())
    del stdout
    del stderr

    return mu_dict


def maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot, f_minarc=1.2):
    """Helper method to grid over B0 at fixed eta2, mu (passed in via pars)

    Procedure:
    1. Improve initial guess for B0 using simple model fit to "average" fwhms
       (midway between fwhms_min, fwhms_max).
    2. From new initial guess for B0, generate grid over many values of B0.
       Store values of B0 and FWHMs from complex model.

    Input:
        pars (lmfit.Parameters): mu, eta2, initial B0 guess, any other parameters needed...
        fwhms_min, fwhms_max must be numpy arrays
    Outputs:
        list of B0 values, list of FWHM value lists (for each B0)
    """

    def rscale(x):  # Grid with rescaled width averaged over FWHMS, r in [0,1]
        return np.mean((x-fwhms_min)/(fwhms_max - fwhms_min))

    def f_rscale(grid_B0):#, rminarc_adapt=True, **kwargs):
        """Rescale model width function, option for adaptive rminarc setting
        kwargs: passed straight to width_cont(...)"""
        pars.add('B0', value=grid_B0)  # Vary B0, other parameters same
        fwhms = width_cont(pars, kevs, snr, rminarc = max(fwhms_max)*f_minarc)
        print '\tModel fwhms = {}'.format(fwhms)
        return rscale(fwhms), fwhms

    # Get FWHMs from initial guess, reinitialize B0 if initial guess is bad
    # WARNING: risk of oscillating forever if function is badly behaved
    print 'Checking initial guess for B0'
    r_init, fwhms_init = f_rscale(pars['B0'].value)
    B0_init = pars['B0'].value
    if all(fwhms_init > fwhms_max):
        pars.add('B0', value=B0_init * 1.1)  # TODO: magic number...
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot)
    elif all(fwhms_init < fwhms_min):
        pars.add('B0', value=B0_init / 1.1)
        return maketab_gridB0(snr, pars, kevs, fwhms_min, fwhms_max, n_tot)

    # Prepare to compute B0 values and FWHMs.
    # Grid evenly over rescaling of FWHMS (rscale)
    n_thin = int(np.around(r_init * n_tot))
    n_wide = int(np.around((1 - r_init) * n_tot))
    dr = 1.0 / n_tot  # y-coord for fill_pts

    print 'Using initial B0 value {} muG'.format(B0_init*1e6)
    print 'Now finding B0 values with FWHMs in range.'
    print 'Require {} values above, {} below initial B0'.format(n_thin, n_wide)

    # Now, actually calculate B0 values.
    if r_init < 1:
        print 'Computing values below initial B0'  # Increasing r towards 1
        B0_wide, r_wide, fwhms_wide = fill_pts(B0_init, fwhms_init, r_init, 1, dr, f_rscale,
                                               dx_max=10e-6)
    else:
        B0_wide, r_wide, fwhms_wide = [B0_init], [r_init], [fwhms_init]
    if r_init > 0:
        print 'Computing values above initial B0'  # Decreasing r towards 0
        B0_thin, r_thin, fwhms_thin = fill_pts(B0_init, fwhms_init, r_init, 0, dr, f_rscale,
                                               dx_max=10e-6)
    else:
        B0_thin, r_thin, fwhms_thin = [B0_init], [r_init], [fwhms_init]
    # Repetitive if/else structure to avoid unnecessary function calls
    # Each call to fill_pts computes the local Jacobian at the same spot (hence
    # the repeated function call -- but, it won't be stored in the output)

    print 'Filling in FWHM range to achieve desired spacing'

    # Combine lists, remove doubly included B0_init/r_init and sort
    B0_all = np.append(B0_thin, B0_wide)
    fwhms_all = np.append(fwhms_thin, fwhms_wide, axis=0)
    r_all = np.append(r_thin, r_wide)

    B0_all, uniqsrt = np.unique(B0_all, return_index=True)
    fwhms_all = fwhms_all[uniqsrt]
    r_all = r_all[uniqsrt]

    # Fill in the gaps
    B0_all, r_all, fwhms_all = mind_the_gaps(B0_all, r_all, fwhms_all, dr, f_rscale)

    return B0_all, fwhms_all  # No longer need r-values


def fill_pts(x0, yaux0, y0, y_final, dy_step, f, dx_max=float('inf'),
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
        likely sorted, or sorted in reverse order
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
               'mu = {:0.3f}; rminarc = {:0.2f}').format(B0*1e6, eta2, mu,
                                                         rminarc)

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


# ======================================================
# Useful tees -- for code-integrated log file generation
# ======================================================
# (in principle we could even parse the logs to recover model FWHMs if
# something in the code goes wrong -- but we don't want to do that...)

class TeeStdout(object):
    """From Python mailing list
    mail.python.org/pipermail/python-list/2007-May/438106.html
    and, stackoverflow.com/a/616686
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
