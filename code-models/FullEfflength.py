"""
A reimplementation of Sean's model code in Python, for clarity?
Also just to make sense of the code

Try to introduce more abstraction

Aaron Tran
2014 July 15
"""

from __future__ import division # Or whatever

import lmfit
import numpy as np
import scipy as sp
import scipy.optimize as opt

# Physical parameters

COMP_RATIO = 4.0  # Rankine-Hugoniot
V0 = 5.0 * 1e8 / compratio  # As measured by Satoru et al., cm/s
R_S = 2.96e19  # tan(943.7") * 2.2 kpc
R_S_ARC = 900.0  # R_S in arcseconds from Green's SNR catalog
ALPHA = 6.0  # Spectral index
S = 2.0*alpha + 1.0

# Constants

C_M = 1.82e18  # For max synchrotron frequency
C_1 = 6.27e18  # For synchrotron emissivity
A = 1.57e-3  # For synchrotron loss time
NUKEV = 2.417989e17  # 1 keV photon frequency

# Options
ICUT = 0  # 1 for cut-off, 0 for no cut-off
IFIT = 1  # 1 for fitting version of code
IPLOT = 0  # 1 for plotting version of code

# Grid parameters
IRMAX = 400  # Resolution on intensity profile
IRADMAX = 100  # Resolution on tabulated electron distribution
IXMAX = 500  # Resolution along line of sight for integration


def main(Bfield, eta2, mu):
    """BLAH"""
    eta = eta2 * (2.0*NUKEV / (C_M * Bfield))**(-(mu-1.0)/2.0)
    # A ton of parameters get set here

    fex = pacholczyk_tab(...)

    for nu in nugraph:
        disttab = tabulate_e_distr(...)  # 2-D np.array in particle energy, r
        jnu = emisx(...)  # either a function or a table over radial pos x

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

# Electron distribution functions
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


# utilities

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
