"""
A reimplementation of Sean's model code in Python, for clarity?
Also just to make sense of the code

Try to introduce more abstraction

Aaron Tran
2014 July 15
"""

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

def emisx(r, nu, B, emis1, radtab, disttab, imu):
    """Numerically integrate 1-particle synchrotron emissivity (emis1) over
    electron distribution to obtain full emissivity j_\nu"""
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
