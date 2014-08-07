"""
Port Sean's model code to Python, aiming for:
1. clarity
2. better speed and resolution

Try to introduce more abstraction

Aaron Tran
2014 July 15
"""

from __future__ import division

import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.integrate as spint

# Constants
C_M = 1.82e18  # For max synchrotron frequency
C_1 = 6.27e18  # For synchrotron emissivity
A = 1.57e-3  # For synchrotron loss time
NUKEV = 2.417989e17  # 1 keV photon frequency


def main():
    pass


# ==============================
# Main function to compute FWHMs
# ==============================

def fefflen(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s, rminarc, icut, irmax,
            iradmax, ixmax):
    """Use numpy arrays for everything (kevs, rminarc)
    
    Input:
        kevs (np.ndarray)
        B0, eta2, mu (float)
        vs, v0, rs, rsarc, s (float)
        rminarc (np.ndarray)
        icut (bool, or 0/1)
        irmax, iradmax, ixmax (int)
            irmax = resolution on intensity profile (400, SN 1006)
            iradmax = resolution on tabulated e- distribution (100, SN 1006)
            ixmax = resolution for line of sight integration (500, SN 1006)
    Output:
        np.ndarray of FWHMs for each entry in kevs
    """
    rmin = 1 - rminarc / rsarc  # array of scaled r, range [0,1]
    # eta = \eta*(E_h)^(1-\mu), per section 4.2 of paper
    eta =  eta2 * (2.*NUKEV/(C_M*B0))**(-(mu-1)/2)

    if icut == 1:
        # Scaling checks out, 1.602177 erg ~= 1 TeV
        # TODO numerical prefactor 8.3 TeV not checked yet
        ecut = ( (8.3*(B0/(100e-6))**(-0.5)*(vs/1e8)*1.6021773)**(2/(1+mu))
                *(1./eta)**(1./(mu+1)) )
    else:
        ecut = 1e40

    widths = []

    for en_nu, rmin_nu in zip(kevs, rmin):
        nu = en_nu*nukev

        # Tabulate e- distr. for fixed nu over e- energy, radial coord
        # TODO build 2-D array efficiently, not point-wise
        xex, fex = read_fglists()
        radtab = np.linspace(rmin_nu, 1.0, num=iradmax)
        disttab = np.empty((len(xex), len(radtab)))
        for i in xrange(len(xex)):
            en = np.sqrt(nu/c1/B0/xex(j))  # e- energy at y = xex(j)
            for irad in xrange(len(radtab)):
                rad = radtab[irad]
                disttab[j,irad] = distr(E, B, ecut, rad, eta, mu, rs, v0, s)

        # rhotab, emistab = ...

        # Compute intensity over rmesh (TODO make adaptive)
        rmesh = np.linspace(rmin_nu, 1.0, num=irmax)
        intensity = np.empty(len(rmesh))
        for i in xrange(rmesh.size):
            r = rmesh[i]

            # Intensity is emissivity integrated over line of sight
            xmesh = np.linspace(0., np.sqrt(1.-r**2), num=ixmax)
            rhomesh = np.sqrt( xmesh**2 + r**2 ) # Even mesh in x, not rho
            emismesh = np.empty(xmesh.size)
            for i in xrange(xmesh.size):
                # Integrate one-particle emissivity over f(x,E) at r=rho
                emismesh[i] = emisx(rhomesh[i], nu, B0, radtab, disttab)
            intensity[i] = spint.simps(emismesh, xmesh)

        dr = (1 - rmin_nu) / (irmax - 1)
        widths.append(fwhm(intensity, dr))


# TODO emisx is called multiple times, and
# we also need xex to tabulate e- distribution
# reading file again and again is super inefficient
# but, keep coding for now and optimize later.

def read_fglists(fname='fglists.dat'):
    """Read frequency-dep. 1-particle emissivity, x and F(x)
    (or, y and G(y) in paper).  Pacholczyk, 1970.
    """
    fgl = np.loadtxt(fname, comments='!')
    if fgl[-1, 0] == 0:  # Truncate Sean's row of zeros at end
        fgl = fgl[:-1,:]
    return fgl[:,0], fgl[:,1]


# =======================================
# Compute emissivity at radial position r
# =======================================

def emisx(r, nu, B, radtab, disttab):
    """Numerically integrate 1-particle synchrotron emissivity (emis1) over
    electron distribution to obtain full emissivity j_\nu

    j_\nu(x) ~ \int_0^\infty y^{-3/2} G(y) f(x, E(y)) dy

    dist = f(x(r), E(y)) as we use r over x, disttab gridded in scaled r

    Inputs:
        r: radial coordinate (scaled on [0,1])
        nu: radiation frequency in Hz
        radtab: table of radial coords, on which electron distr is sampled
                (set by iradmax; compare irmax for intensity tables)
        disttab: electron distribution function, gridded as:
                 xex(j) (particle energy / radiation frequency), radial coord
    Output: j_\nu(r), r = radial dist, j_\nu a scalar (float)
    """
    xex, fex = read_fglists()

    # TODO use scipy.integrate for speed up

    # Initialize integral (trapezoidal sum)
    dist =  interp(r, radtab, disttab[0,:])
    x_old = xex[0]
    intgd_old = fex[0] * dist / (xex[0]**1.5)
    trpsum = 0.
    
    # Integrate over y-values (e- energies) in Pacholczyk
    for i in xrange(1, len(xex)):
        dist = interp(r, radtab, disttab[i,:])

        x = xex[i]
        intgd = fex[i] * dist / (x**1.5)
        trpsum += (intgd + intgd_old) * (x - xold)

        x_old = x
        intg_old = intg

    return trpsum * np.sqrt(nu * B)  # Constant prefactors skipped


# ===============================
# Electron distribution functions
# ===============================

# Not functional -- subroutines need access to constants

def distr(E, B, Ecut, r, eta, mu, rs, v, s):
    """Electron distribution, from Rettig & Pohl and Lerche & Schlickeiser

    Input:
        E = particle energy, B = magnetic field
        ecut = cut-off energy of injected electrons
        r = radial position (scaled to [0,1])
        eta = \eta*(E_h)^(1-\mu), from the paper
        mu = diffusion coefficient exponent
        rs, v, s = shock radius, plasma velocity, injected e- spectral index
    """

    D0 = eta*2.083e19/B  # eta (E_h)^(1-mu) * C_d / B0
    a = 1.57e-3  # This is b in the paper
                 # Synchrotron constant ... TODO declare at module level?
    k0 = 1.  # Normalization, Q_0
    b0 = a * B**2  # b*B^2 in paper
    alpha = 1. # TODO Not sure what this is, scaling factor?

    tau = 1. / (b0 * E)  # Synchrotron cooling time
    lad = v * tau  # Advective lengthscale
    ldiff = np.sqrt(D0*E**mu*tau)  # Diffusive lengthscale

    x = (rs - r*rs)/alpha  # Convert scaled r to x in cm

    # lad/ldiff ~ sqrt(v*lad) / dsqrt(D) ~ peclet number!
    if lad/ldiff > 30:
        return distr_adv(E, B, Ecut, r, eta_c, mu, rs, v, s)

    if mu > 1:
        return distr_mgt1(E, B, Ecut, r, eta_c, mu, rs, v, s)
    elif mu < 1:
        return distr_mlt1(E, B, Ecut, r, eta_c, mu, rs, v, s)
    else:
        return distr_rpohl(E, B, Ecut, r, eta_c, mu, rs, v, s)

# TODO exponentiation (**) works elementiwse on numpy arrays 

def distr_adv(E, B, Ecut, r, eta_c, mu, rs, v, s):
    """Case: l_ad / l_diff >> 1"""
    efinv = a/v * B**2 * (rs - r*rs)
    ef = 1. / efinv
    en0 = E / (1.0 - E/ef)
    argexp = en0/ecut
    # argexp = E / ecut / (1 - a*B**2*E * (rs-r*rs) / v)
    Xi = 1 - E/ef  # Only used for if-else structure below
    if Xi > 0:
        return 1./en0**s * (en0/E)**2 / np.exp(argexp) * 8e-9
    else
        return 0.


def distr_mlt1(E, B, Ecut, r, eta_c, mu, rs, v, s):
    """Case: mu < 1

    Integrate over t, where n = 1 - t^2
    """
    t = np.linspace(0., 1., num=1000)  # TODO arbitrary
    n = 1 - t**2

    argexp = (n**(-1/(1-mu)) * E/Ecut +
                  (lad/alpha * (1 - n**(1/(1-mu)) ) - x)**2 /
                  (4*ldiff**2*(1-n)/alpha) * (1-mu) )
    argexp[n==1.] = 1e35  # Prevent blowup
    integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)

    integral = spint.trapz(integrand, t)  # TODO can do better?

    return integral * ( k0 /2 /np.sqrt(pi*alpha*b0*D0*(1-mu)) *
                        E**(-1*(mu/2+1/2+s)) )


def distr_rpohl(E, B, Ecut, r, eta_c, mu, rs, v, s):
    """Case: mu = 1
    Integrate over t in [0,1] (n=e^(t^2) in [1,e]), then
    integrate over y in [0,1] (n=mess in [e, infty)).
    """
    t = np.linspace(0., 1., num=1000)  # TODO arbitrary
    n = np.exp(t**2)

    argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
                          (4*ldiff**2/alpha*np.log(n)) )
    argexp[np==1.] = 1e35  # Prevent blowup
    integrand = 2*np.exp((1-s)*t**2) / np.exp(argexp)
    integral_left = spint.trapz(integrand, t)  # TODO can do better?

    y = np.linspace(0., 1., num=100)  # TODO arbitrary
    q = 2.
    nmin = np.exp(1.)
    n = y / (1 - y**2)**q + nmin

    argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
                          (4*ldiff**2/alpha*np.log(n)) )  # same
    integrand = ( n**(-1*s) / np.sqrt(np.log(n)) / np.exp(argexp) *
                  ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
    integral_right = spint.trapz(integrand, y)  # TODO can do better?

    integral = integral_left + integral_right
    return integral * k0 /2 /np.sqrt(pi*alpha*b0*D0) * E**(-1*(1+s))


def distr_mgt1(E, B, Ecut, r, eta_c, mu, rs, v, s):
    """Case: mu > 1
    Integrate over t in [0, 1] (n = 1+t^2 in [1, 2]), then
    integrate over y in [0, 1] (n = mess in [2, infty)).
    """
    t = np.linspace(0., 1., num=1000)  # TODO arbitrary
    n = 1 + t**2
    argexp = n**(-1/(1-mu)) * E/Ecut + ( (lad/alpha*(1-n**(1/(1-mu)))-x)**2 /
                                         (4*ldiff**2*(1-n)/alpha) * (1-mu) )
    argexp[n==1.] = 1e35
    integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)
    integral_left = spint.trapz(integrand, t)  # TODO can do better?

    y = np.linspace(0., 1., num=100)  # TODO arbitrary
    q = 2.
    nmin = 2.
    n = y / (1 - y**2)**q + nmin
    argexp = ( n**(-1/(1-mu)) * E/Ecut +
               (lad/alpha*(1-n**(1/(1-mu))) - x)**2 /
               (4*ldiff**2*(1-n)/alpha) * (1-mu) )
    integrand = ( n**((s+mu-2)/(1-mu)) / np.sqrt(n-1) / np.exp(argexp) *
                  ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
    integral_right = spint.trapz(integrand, t)  # TODO can do better?

    integral = integral_left + integral_right
    return integral * ( k0 /2 /np.sqrt(pi*alpha*b0*D0*(mu-1)) *
                        E**(-1*(mu/2+1/2+s)) )


# ============================
# Find FWHM in intensity graph
# ============================

def fwhm(rmesh, intensity, dr, nu):
    """Compute fwhm. nu is for error printing"""

    # Compute half max and find crossings
    half_max = 0.5 * np.amax(intensity)
    cross = np.diff(np.sign(intensity - half_max))
    inds_rmin = np.where(cross > 0)[0]  # Left (neg to pos)
    inds_rmax = np.where(cross < 0)[0]  # Right (pos to neg)

    if inds_rmin.size == 0:  # No left crossing found
        print 'Box length error (rmin) at {} keV'.format(nu/NUKEV)
        width = 1.  # Maximum (scaled) width
        return width
    else:
        i_rmin = inds_rmin[-1]  # Crossing closest to peak

    if inds_rmax.size == 0:  # No right crossings found
        print 'Box length error (rmax) at {} keV'.format(nu/NUKEV)
        print 'Rim falloff towards shock is weird?'
        i_rmax = -1  # Farthest right index
    else:
        i_rmax = inds_rmax[0] + 1  # Crossing closest to peak

    width = rmesh(i_rmax) - rmesh(i_rmin)

    if width / abs(dr) < 20:  # Arbitrary resolution constrain
        print 'Resolution error at {} keV'.format(nu/NUKEV)
        width = 0.

    return width


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
