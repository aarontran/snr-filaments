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

import lmfit
import matplotlib.pyplot as plt
import numpy as np
from numpy import f2py
import scipy as sp
from scipy import optimize

import snr_catalog as snrcat

# Compile and import Fortran code -- not sure where to do this.
# Could throw this into another short script? remember to call
# fm.readfglists() before anything else though.
with open('FullEfflength_mod.f', 'r') as f:
    fsource = ''.join(list(f))
np.f2py.compile(fsource, modulename='fullmodel', verbose=1)
import fullmodel as fm
fm.readfglists()  # Only do this once


def main():
    """Just testing, for now"""

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
    p.add('B0', value=1e-4, min=1e-6, max=1e-2)
    #p.add('eta2', value=1, vary=False)
    p.add('eta2', value=1, min=0, max=1e4)
    p.add('mu', value=1, vary=False)


    # Fit complex model
    res = lmfit.minimize(objectify(width_cont), p,
                   args=(data, eps, kevs, snr),
                   kwargs={icut:0},  # Numerical code settings
                   method='leastsq',
                   epsfcn=1e-6)  # this epsfcn is IMPORTANT!!!

    # Fit simple model
    res = lmfit.minimize(objectify(width_cont), p,
                   args=(data, eps, kevs, snr),
                   method='leastsq')

    # Print fit information
    lmfit.printfuncs.report_fit(res.params)
    ci = lmfit.conf_interval(res)  # Only if more than 2 variables?
    lmfit.printfuncs.report_ci(ci)

    # Plot fit with data
    plot_data()
    en_plt = np.linspace(0.6, 2.1, 10)
    plt.errorbar(kevs, data, eps, fmt='bo')
    plt.plot(en_plt, width_cont(p, en_plt), '-k.')
    plt.xlabel('Observation energy (keV)')
    plt.ylabel('FWHM (arcsec)')
    plt.show()

    # Maybe good to spit out a configuration file here, to log fit inputs and
    # outputs.

    return res


# ==========================================
# Functions to execute various types of fits
# ==========================================

def fit_from_tab():
    """Fit data from precomputed table of values..."""
    pass


def make_tab():
    """SNR dependent tabulation of FWHM values for fitting / explore parameter space"""
    pass


def manual_fit_a_la_sean(snr):
    """Interactively prompts user for values of B0, eta2, mu
    Prints FWHMs, residuals, and chi-squared values in the process
    To replicate work in SN 1006, use only 2 values, and computes a FWHM and an
    exponent mE (both at 2keV) for manual fitting.
    """
    pass


def manual_fit(snr):
    """Interactively prompts user for values of B0, eta2, mu
    Prints FWHMs, residuals, and chi-squared values in the process
    """
    pass


# =============================
# Model functions, objectify(f)
# =============================

def width_cont(params, kevs, snr, rminarc=60, icut=1, irmax=400, iradmax=100, ixmax=500):
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

    print ('Function call with B0 = {:0.3f} muG; eta2 = {:0.3f}; '
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
    return lambda pars, data, eps_data, *args, **kwargs:
        (f(pars, *args, **kwargs) - data)/eps_data


if __name__ == '__main__':
    main()
