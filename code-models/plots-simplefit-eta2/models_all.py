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
# 2. tab_eta2(...) - 
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

    This code traces out a swath in the eta2-B0 plane (approx. 1-D curve)
    but does not grid over mu

    Need to configure.
    
    May need to allow to pick eta2 gridding parameters (range of values, number
    of data points).

    Finally, think about how to store the data -- we need to store the 3-tuple
    of eta2, B0, mu, along with model FWHMs (don't recompute each time you try
    to perform a fit).

    Output values of chisqr_vals are SPECIFIC to this particular set of
    FWHMs...

    I have a single B0 value for each pair (eta2, mu)
    but we need a few values of B0, to span possible FWHM range...
    Time to implement B0 rangefinding procedure...

    so for each mu (only one in this method),
           for each eta2, return B0 values + model FWHMs for each B0

    use the input FWHMs as a guideline/baseline to get simple model fit values.
    Regrettably, the tabulation will be FWHM input dependent...
    But, in the end, the result should not depend strongly on the inputs.
        YOU NEED TO VERIFY THIS -- generate TWO TABLES, with slightly different
        measured FWHM inputs (say, two different filaments).
        Then, check that output results from full procedure are consistent
        within error whatever...

        ANDDDDDDD YOU need to think about how to get fricking error from this
        procedure.  I think once I have the best fit B0, eta2, mu then I can
        deal with error and that shouldn't be so bad.

    Friday -- bring tables of results to brian (cc rob)
        bring poster sketch / draft to discuss things to include, briefly...
    """

    # For set of eta2 values, find best fit B0 values from simple model
    # B0 values from simple model are reasonably close to complex model
    # but complex model requires smaller B0 (at fixed eta2)

    eta2_vals = np.logspace(-2, 3, 10, base=10) # If using full code
    #eta2_vals = np.logspace(-3, 4, 200, base=10) # plotting (simple only)

    # Set-up for loop
    p = lmfit.Parameters()
    p.add('mu', value=mu, vary=False)

    b0_vals = []
    chisqr_vals = []
    prev_b0 = 1e-6

    full_model_fwhms = []  # List of lists/arrays

    np.set_printoptions(precision=2)

    for eta2 in eta2_vals:
        p.add('eta2', value=eta2, vary=False)
        p.add('B0', value=150e-6, min = 1e-6, max=1e-2)  # Just initial guess

        print '\nGrid: eta2 = {:.2f}, mu = {:.2f}'.format(eta2, mu)
        print '--------------------------------'

        # Run simple model fit
        res = lmfit.minimize(objectify(width_dump), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq')
        b0_simp = p['B0'].value
        chi2_simp = res.chisqr
        fwhms_simp = width_dump(p, kevs, snr)

        # Here I need to generate multiple values of B0, which work.
        # Need to save the values that I do manage to tabulate...
        # (no sense in wasting generated intermediate values)

        # here I'm expressly NOT executing a best fit -- a fit may not
        # hit all the desired values / sample them well.

        # Run complex model fit
        # Use previous complex fit B0, and current simple fit B0, to help
        # constrain current fit B0
        if prev_b0 == 1e-6:
            p.add('B0', value=p['B0'].value, min=prev_b0*0.8, max=p['B0'].value*1.5)
        else:
            p.add('B0', value=prev_b0, min=prev_b0*0.8, max=p['B0'].value*1.5)

        print 'Simple fit B0 = {:.2f} muG'.format(b0_simp*1e6)
        print 'Starting complex model fit (eta2 fixed, B0 free)'
        res = lmfit.minimize(objectify(width_cont), p,
                             args=(data, eps, kevs, snr),
                             kws={'icut':1},
                             method='leastsq',
                             epsfcn=1e-1)  # changed from 1e-5, due to B0 bounds
        b0_comp = p['B0'].value
        chi2_comp = res.chisqr
        fwhms_comp = width_cont(p, kevs, snr, verbose=False)
        prev_b0 = b0_comp # To help constrain next fit
        # since next eta2 is bigger, next B0 must increase as well

        print 'measured:    {}'.format(data)
        print ('simple fit: {}; '.format(fwhms_simp)
                + 'B0 = {:.2f} muG, chisqr = {:.2f}'.format(b0_simp*1e6, chi2_simp))
        print ('full fit:   {}; '.format(fwhms_comp)
               + 'B0 = {:.2f} muG, chisqr = {:.2f}'.format(b0_comp*1e6, chi2_comp))

        #b0_vals.append(b0_simp)
        #chisqr_vals.append(chi2_simp)
        b0_vals.append(b0_comp)
        chisqr_vals.append(chi2_comp)
        full_model_fwhms.append(fwhms_comp)


    # Additional fit to find best fit eta2, B0 in tandem
    # but this should be run separately -- doing a slightly different job.
    # This is important for reproducing Table 7 of Ressler!
    # And, this should agree with the output plots from this tabulation process
    p.add('B0', value=150e-6, min=1e-6, max=1e-2)
    p.add('eta2', value=1, min=0.01, max=1000)
    res = lmfit.minimize(objectify(width_dump), p,
                             args=(data, eps, kevs, snr),
                             method='leastsq')
    print 'Best fit values with both eta2, B0 free'
    print 'eta2 = {:.2f}, B0 = {:.2f}, mu = {:.2f}'.format(p['eta2'].value,
                                                     p['B0'].value * 1e6,
                                                     p['mu'].value)
    print '---------------------------------------'

    print 'measured: {}'.format(data)
    print ('simple fit: {}; '.format(width_dump(p, kevs, snr))
            + 'chisqr = {:.2f}'.format(res.chisqr))

    return eta2_vals, b0_vals, chisqr_vals, fwhm_vals


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

def width_cont(params, kevs, snr, verbose=True, rminarc=60, icut=1,
               irmax=400, iradmax=100, ixmax=500):
    """Width function, wraps numerical computation in fm.fullefflengthsub
    which solves continuous loss model (equation 12)

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
