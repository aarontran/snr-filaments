"""
A reimplementation of Sean's model code in Python, for clarity?
Also just to make sense of the code

Try to introduce more abstraction

Aaron Tran
2014 July 15
"""

from __future__ import division # Or whatever

import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.integrate as spint

C_M = 1.82e18  # For max synchrotron frequency
C_1 = 6.27e18  # For synchrotron emissivity
A = 1.57e-3  # For synchrotron loss time
NUKEV = 2.417989e17  # 1 keV photon frequency


def main():
    pass


# ==============================
# Main function to compute FWHMs
# ==============================

def fefflen(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s, rminarc,
            icut, irmax, iradmax, ixmax):
    """LALALA. Use numpy arrays for everything"""

    # Calculate derived parameters
    rmin = 1 - rminarc / rsarc  # vector
    eta =  eta2 * (2.*NUKEV/(C_M*B0))**(-(mu-1)/2)

    if icut == 1:
        ecut = (8.3*(B0/(100e-6))**(-0.5)*(vs/1e8)
                *1.6021773)**(2./(1.+mu)) * (1./eta)**(1./(mu+1))
    else:
        ecut = 1e40

    # Start integrating
    for en in kevs:
        nu = en*nukev

        radtab, disttab = get_disttab(...)
        # rhotab, emistab = get_emisx(radtab, disttab, ...)
        # NOTE Think about sampling here... how fine is needed?

        # Integrate over line of sight
        # Worry about "smart" / adaptive meshing later...
        for r in rmesh:

            rhomesh = stuff
            emisrho = emisx(rhomesh)  # Integrate over particle energy f(x,E)
                                      # giving emissivity at each point
            # emisrho = map(lambda rho: interp(rho, rhotab, emistab), rhomesh)
            intensity = sp.integrate.quad(...) # Integrate over x

        widths = find_fwhm(intensities)

# Functions to compute integrals, tabulate stuff

def emisx(r, nu, B, radtab, disttab):
    """Numerically integrate 1-particle synchrotron emissivity (emis1) over
    electron distribution to obtain full emissivity j_\nu
    Inputs:
        r: radial coordinate (TODO: scaled, arcsec, ??)
        nu: radiation frequency in Hz
        radtab: table of radial coords, on which electron distr is sampled
                (set by iradmax; compare irmax for intensity tables)
        disttab: electron distribution function, gridded as:
                 xex(j) (particle energy / radiation frequency), radial coord
    Output: j_\nu(r), r = radial dist, j_\nu a scalar (float)
    """
    pass


# ===============================
# Electron distribution functions
# ===============================

# throw this into a CLASS or something -- so I can hold CONSTANTS
# at a class level or whatever...

def distr(E, B, E_cut, r, eta_c, mu, rs, v, s):
    """Electron distribution, from Rettig & Pohl and Lerche & Schlickeiser"""

    D0 = eta*2.083d19/B  # CHECK THIS eta (E_h)^(1-mu) * C_d / B0
    a = 1.57e-3  # Synchrotron constant ... declare at module level
    tau = 1. / (a * B**2 * E)
    lad = v * tau
    ldiff = dsqrt(D0*E**mu*tau)

    # lad/ldiff ~ sqrt(v*lad) / dsqrt(D) ~ peclet number!

    if lad/ldiff > 30:
        return distr_adv()

    if mu > 1:
        return distr_mgt1(E, B, E_cut, r, eta_c, mu, rs, v, s)
        # I think this is the same as below...
    elif mu < 1:
        return distr_mlt1(E, B, E_cut, r, eta_c, mu, rs, v, s)
    else:
        return distr_rpohl(E, B, E_cut, r, eta_c, mu, rs, v, s)


def distr_adv(E, B, E_cut, r, eta_c, mu, rs, v, s):
    """If advective length much greater than diffusive lengthscale"""
    pass


def distr_mgt1(E, B, E_cut, r, eta_c, mu, rs, v, s):
    """Case: mu > 1"""
    pass


def distr_mlt1(E, B, E_cut, r, eta_c, mu, rs, v, s):
    """Case: mu < 1"""
    pass


def distr_rpohl(E, B, E_cut, r, eta_c, mu, rs, v, s):
    """Case: mu = 1"""
    pass


# =========
# Utilities
# =========

# ALTERNATELY, just use np.interp...

def interp(x, xdat, ydat):
    """Linear interpolation of data at x; xdat MUST be sorted (for speed)"""
    ind = np.searchsorted(xdat, x)  # xdat[ind-1] < x < xdat[ind]
    return interp_ind(x, ind, xdat, ydat)

def interp_ind(x, ind, xdat, ydat):
    """Linear interpolation of data at x, xdat[ind-1] < x < xdat[ind]
    (premature) optimization; avoid repeated binary search"""
    slope = (ydat[ind] - ydat[ind-1]) / (xdat[ind] - xdat[ind-1])
    return ydat[ind-1] + slope * (x - xdat[ind-1])


if __name__ == '__main__':
    main()
