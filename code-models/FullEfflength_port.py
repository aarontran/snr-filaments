"""
Port Sean's model code to Python, aiming for:
1. clarity
2. better speed and resolution

Aaron Tran
2014 August 6-8
"""

from __future__ import division

from datetime import datetime
import numpy as np
import scipy as sp
import scipy.optimize as spopt
import scipy.integrate as spint

import fullmodel
fullmodel.readfglists()  # !! important

# Constants
C_M = 1.82e18  # For max synchrotron frequency
C_1 = 6.27e18  # For synchrotron emissivity
A = 1.57e-3  # For synchrotron loss time
NUKEV = 2.417989e17  # 1 keV photon frequency


def main():
    #B0 = raw_input('Enter B0 (Gauss): ')
    #eta2 = raw_input('Enter eta2: ')
    #mu = raw_input('Enter mu: ')
    B0 = 150e-6
    eta2 = 1
    mu = 1.5

    kevs = np.array([0.7, 1., 2.])
    rminarc = 20. * np.ones(3)

    fwhms = fefflen(kevs, B0, eta2, mu, 5e8, 5e8/4., 2.96e19, 900.,
                    2.2, rminarc, True, 100, 200, 1000)
    
    for en, fwhm in zip(kevs, fwhms):
        print '{:0.2f} keV: {:0.5f}'.format(en, fwhm)

    return fwhms

# ==============================
# Main function to compute FWHMs
# ==============================

# @profile
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

    # Various derived constants
    rmin = 1 - rminarc / rsarc  # array of scaled r, range [0,1]
    # eta = \eta*(E_h)^(1-\mu), per section 4.2 of paper
    eta =  eta2 * (2.*NUKEV/(C_M*B0))**(-(mu-1)/2)
    if icut:
        Ecut = ( (8.3*(B0/(100e-6))**(-0.5)*(vs/1e8)*1.6021773)**(2/(1+mu))
                *(1./eta)**(1./(mu+1)) )  # TODO check prefactor 8.3 TeV
    else:
        Ecut = 1e40

    # 1-particle synchrotron emissivity
    xex, fex = read_fglists()

    widths = []

    for en_nu, rmin_nu in zip(kevs, rmin):
        nu = en_nu*NUKEV

        # Tabulate e- distr. for fixed nu over e- energy, radial coord
        # Be careful, Fortran code hardcoded for iradmax < 1000
        # and Pacholczyk table < 100 entries
        #print 'e- tabulation start: {}'.format(datetime.now())
        radtab = np.linspace(rmin_nu, 1.0, num=iradmax)
        disttab = fullmodel.distr(np.zeros(1000), iradmax, nu, rmin_nu,
                                  B0, Ecut, eta, mu, rs, v0, s, C_1)
        disttab = disttab[:len(xex), :iradmax]  # Slice off zeros
        #print 'e- tabulation done: {}'.format(datetime.now())

        # Tabulate emissivity over radial coord
        #print 'emis tab start: {}'.format(datetime.now())
        rhotab = np.linspace(rmin_nu, 1., num=10000)  # TODO set as irhomax?
        emistab = emisx(rhotab, nu, B0, radtab, disttab, xex, fex)
        #print 'emis tab done: {}'.format(datetime.now())

        # Compute intensity over rmesh
        #print 'intensity start: {}'.format(datetime.now())
        rmesh = np.linspace(rmin_nu, 1.0, num=irmax, endpoint=False)
        xmaxmesh = np.sqrt(1.-rmesh**2)  # Max x on line of sight at each r
        # Grid along line of sight for each r; shape is (irmax, ixmax)
        # [[0., ..., xmax_0], [0., ..., xmax_1], ... [0., ..., xmax_{irmax-1}]]
        xintv = np.linspace(0., 1., num=ixmax)  # Line of sight (x) sampling
        xgrid = xintv * xmaxmesh.reshape(-1,1)  # stackoverflow.com/a/16887425
        # Compute rho, emissivity at each point in grid
        rgrid = np.tile(rmesh, (len(xintv), 1)).T  # r = constant on sightlines
        rhogrid = np.sqrt(xgrid**2 + rgrid**2)
        emisgrid = np.interp(rhogrid, rhotab, emistab)
        # then integrate emis. on lines of sight
        intensity = spint.simps(emisgrid, xgrid)
        #print 'intensity done: {}'.format(datetime.now())

        def get_intensity(r):
            """Non-vectorized intensity computation, for FWHM finding"""
            xmesh = np.linspace(0., np.sqrt(1.-r**2), num=ixmax)
            rhomesh = np.sqrt(xmesh**2 + r**2)
            emismesh = np.interp(rhomesh, rhotab, emistab)
            return spint.simps(emismesh, xmesh)

        #print 'find fwhm start: {}'.format(datetime.now())
        dr = (1 - rmin_nu) / (irmax - 1) # Match rmesh
        #w = fwhm(rmesh, intensity, dr, nu)
        w = fwhm2(rmesh, intensity, get_intensity, dr, nu)
        widths.append(w)
        #print 'find fwhm done: {}\n'.format(datetime.now())

    return np.array(widths) * rsarc


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

def emisx(r, nu, B, radtab, disttab, xex, fex):
    """Numerically integrate 1-particle synchrotron emissivity (emis1) over
    electron distribution to obtain full emissivity j_\nu

    j_\nu(x) ~ \int_0^\infty y^{-3/2} G(y) f(x, E(y)) dy

    dist = f(x(r), E(y)) as we use r over x, disttab gridded in scaled r

    Inputs:
        r: radial coordinate, array (scaled on [0,1])
        nu: radiation frequency in Hz
        radtab: table of radial coords, on which electron distr is sampled
                (set by iradmax; compare irmax for intensity tables)
        disttab: electron distribution function, gridded as:
                 xex(j) (particle energy / radiation frequency), radial coord
    Output: j_\nu(r), r = radial dist, j_\nu a scalar (float)
    """
    # Integrate over y-values (e- energies) in Pacholczyk
    ind = np.searchsorted(radtab, r)
    slp = (disttab[:,ind] - disttab[:,ind-1]) / (radtab[ind] - radtab[ind-1])
    dist = disttab[:,ind-1] + slp * (r - radtab[ind-1])
    # Would be nice to use np.interp, but not sure if possible

    xex = np.tile(xex, (len(r),1)).T
    fex = np.tile(fex, (len(r),1)).T

    integrand = fex * dist / (xex**1.5)

    # Skip constant prefactors
    return spint.trapz(integrand, xex, axis=0) * np.sqrt(nu*B)


# ===============================
# Electron distribution functions
# ===============================

# NOTE not currently in use!  Calling Fortran code instead

def distr(Ec, rc, B, Ecut, eta, mu, rs, v, s):
    """Electron distribution, from Rettig & Pohl and Lerche & Schlickeiser

    ASSUMED: E, r are 1-D arrays
    spint.trapz integrates along last axis by default, conveniently

    Input:
        E = particle energy, r = radial position (scaled to [0,1])
        B = magnetic field
        Ecut = cut-off energy of injected electrons
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

    def get_nums(E, r):
        """Compute derived parameters over 3-D arrays from np.meshgrid"""
        tau = 1. / (b0 * E)  # Synchrotron cooling time
        lad = v * tau  # Advective lengthscale
        ldiff = np.sqrt(D0*E**mu*tau)  # Diffusive lengthscale
        x = (rs - r*rs)/alpha  # Convert scaled r to x in cm
        return tau, lad, ldiff, x

    # Functions for each of 4 cases

    # TODO this one is not yet functional
    def distr_adv(Ec, rc, B, Ecut, eta, mu, rs, v, s):
        """Case: l_ad / l_diff >> 1"""
        efinv = a/v * B**2 * (rs - r*rs)
        ef = 1. / efinv
        en0 = E / (1.0 - E/ef)
        argexp = en0/Ecut
        # argexp = E / Ecut / (1 - a*B**2*E * (rs-r*rs) / v)
        Xi = 1 - E/ef  # Only used for if-else structure below
        if Xi > 0:
            return 1./en0**s * (en0/E)**2 / np.exp(argexp) * 8e-9
        else:
            return 0.


    def distr_mlt1(Ec, rc, B, Ecut, eta, mu, rs, v, s):
        """Case: mu < 1
        Integrate over t, where n = 1 - t^2
        """
        # Build 3-D arrays
        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
        n = 1 - t**2
        tau, lad, ldiff, x = get_nums(E, r)

        argexp = (n**(-1/(1-mu)) * E/Ecut +
                      (lad/alpha * (1 - n**(1/(1-mu)) ) - x)**2 /
                      (4*ldiff**2*(1-n)/alpha) * (1-mu) )
        argexp[n==1.] = 1e35  # Prevent blowup
        integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)

        integral = spint.trapz(integrand, t)  # TODO can do better?
        E_s = E[:,:,0]  # Slice as integral is a 2-D array
        return integral * ( k0 /2 /np.sqrt(np.pi*alpha*b0*D0*(1-mu)) *
                            E_s**(-1*(mu/2+1/2+s)) )


    def distr_rpohl(Ec, rc, B, Ecut, eta, mu, rs, v, s):
        """Case: mu = 1
        Integrate over t in [0,1] (n=e^(t^2) in [1,e]), then
        integrate over y in [0,1] (n=mess in [e, infty)).
        """
        # First integral, build 3-D arrays
        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
        n = np.exp(t**2)
        tau, lad, ldiff, x = get_nums(E, r)

        argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
                              (4*ldiff**2/alpha*np.log(n)) )
        argexp[n==1.] = 1e35  # Prevent blowup
        integrand = 2*np.exp((1-s)*t**2) / np.exp(argexp)
        integral_left = spint.trapz(integrand, t)  # TODO can do better?

        # Second integral, build 3-D arrays
        yc = np.linspace(0., 1., num=100, endpoint=False)  # TODO arbitrary
        E, r, y = np.meshgrid(Ec, rc, yc, indexing='ij')
        q = 2.  # q, nmin just for change of variables
        nmin = np.exp(1.)
        n = y / (1 - y**2)**q + nmin
        tau, lad, ldiff, x = get_nums(E, r)

        argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
                              (4*ldiff**2/alpha*np.log(n)) )  # same
        integrand = ( n**(-1*s) / np.sqrt(np.log(n)) / np.exp(argexp) *
                      ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
        integral_right = spint.trapz(integrand, y)  # TODO can do better?

        # Final integral
        integral = integral_left + integral_right
        E_s = E[:,:,0]  # Slice as integral is 2-D array
        return integral * k0 /2 /np.sqrt(np.pi*alpha*b0*D0) * E_s**(-1*(1+s))


    def distr_mgt1(Ec, rc, B, Ecut, eta, mu, rs, v, s):
        """Case: mu > 1
        Integrate over t in [0, 1] (n = 1+t^2 in [1, 2]), then
        integrate over y in [0, 1] (n = mess in [2, infty)).
        """
        # First integral, build 3-D arrays
        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
        n = 1 + t**2
        tau, lad, ldiff, x = get_nums(E, r)

        argexp = n**(-1/(1-mu)) * E/Ecut + ( (lad/alpha*(1-n**(1/(1-mu)))-x)**2 /
                                             (4*ldiff**2*(1-n)/alpha) * (1-mu) )
        argexp[n==1.] = 1e35
        integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)
        integral_left = spint.trapz(integrand, t)  # TODO can do better?

        # Second integral, build 3-D arrays
        yc = np.linspace(0., 1., num=100, endpoint=False)  # TODO arbitrary
        E, r, y = np.meshgrid(Ec, rc, yc, indexing='ij')
        q = 2.
        nmin = 2.
        n = y / (1 - y**2)**q + nmin
        tau, lad, ldiff, x = get_nums(E, r)

        argexp = ( n**(-1/(1-mu)) * E/Ecut +
                   (lad/alpha*(1-n**(1/(1-mu))) - x)**2 /
                   (4*ldiff**2*(1-n)/alpha) * (1-mu) )
        integrand = ( n**((s+mu-2)/(1-mu)) / np.sqrt(n-1) / np.exp(argexp) *
                      ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
        integral_right = spint.trapz(integrand, y)  # TODO can do better?

        integral = integral_left + integral_right
        E_s = E[:,:,0]  # Slice as integral is a 2-D array
        return integral * ( k0 /2 /np.sqrt(np.pi*alpha*b0*D0*(mu-1)) *
                            E_s**(-1*(mu/2+1/2+s)) )

    # Need to address edge case by vectorizing... or something.
    # Construct 2-D grid of E, r

    # lad/ldiff ~ sqrt(v*lad) / dsqrt(D) ~ peclet number!
    #if lad/ldiff > 30:
    #    return distr_adv(E, r, B, Ecut, eta, mu, rs, v, s)

    if mu > 1:
        return distr_mgt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
    elif mu < 1:
        return distr_mlt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
    else:
        return distr_rpohl(Ec, rc, B, Ecut, eta, mu, rs, v, s)


# =========================
# Attempt to invoke f2py...
# =========================

def distr2(Ec, rc, B, Ecut, eta, mu, rs, v, s):
    """Electron distribution, from Rettig & Pohl and Lerche & Schlickeiser

    ASSUMED: E, r are 1-D arrays
    spint.trapz integrates along last axis by default, conveniently

    Input:
        E = particle energy, r = radial position (scaled to [0,1])
        B = magnetic field
        Ecut = cut-off energy of injected electrons
        eta = \eta*(E_h)^(1-\mu), from the paper
        mu = diffusion coefficient exponent
        rs, v, s = shock radius, plasma velocity, injected e- spectral index
    """
    disttab = np.empty(Ec.size, rc.size)

    if mu > 1:
        return distr_mgt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
    elif mu < 1:
        return distr_mlt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
    else:
        return distr_rpohl(Ec, rc, B, Ecut, eta, mu, rs, v, s)


# ============================
# Find FWHM in intensity graph
# ============================

def fwhm(rmesh, intensity, dr, nu):
    """Compute fwhm. nu is for error printing."""

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
        i_rmin = inds_rmin[-1]+1  # Crossing closest to peak

    if inds_rmax.size == 0:  # No right crossings found
        print 'Box length error (rmax) at {} keV'.format(nu/NUKEV)
        print 'Rim falloff towards shock is weird?'
        i_rmax = -1  # Farthest right index
    else:
        i_rmax = inds_rmax[0]  # Crossing closest to peak

    width = rmesh[i_rmax] - rmesh[i_rmin]

    if width / abs(dr) < 20:  # Arbitrary resolution constrain
        print 'Resolution error at {} keV'.format(nu/NUKEV)
        width = 0.

    return width


def fwhm2(rmesh, intensity, f_int, dr, nu):
    """Compute fwhm but better (test). nu is for error printing."""

    eps2 = np.finfo(float).eps * 2  # Brentq tolerance

    # Compute half max from grid initial guess
    idxmax = np.argmax(intensity)
    res = spopt.minimize_scalar(lambda x: -1*f_int(x), method='bounded',
                                bounds=(rmesh[idxmax-1], rmesh[idxmax+1]),
                                options={'xatol':eps2})
    halfpk = 0.5 * f_int(res.x)

    # Find where FWHMs should be 
    cross = np.diff(np.sign(intensity - halfpk))
    inds_rmin = np.where(cross > 0)[0]  # Left (neg to pos)
    inds_rmax = np.where(cross < 0)[0]  # Right (pos to neg)

    # Throw errors if we can't identify FWHM
    if inds_rmin.size == 0:  # No left crossing found
        print 'Box length error (rmin) at {} keV'.format(nu/NUKEV)
        return 1. # Maximum (scaled) width
    else:
        i_rmin = inds_rmin[-1]+1  # Crossing closest to peak

    if inds_rmax.size == 0:  # No right crossings found
        print 'Box length error (rmax) at {} keV'.format(nu/NUKEV)
        print 'Rim falloff towards shock is weird?'
        i_rmax = -1  # Farthest right index
        return rmesh[i_rmax] - rmesh[i_rmin]
    else:
        i_rmax = inds_rmax[0]  # Crossing closest to peak

    # Finally, nail down the fwhm location
    def f_thrsh(r):  # for brentq rootfinding
        return f_int(r) - halfpk
    rmin = spopt.brentq(f_thrsh, rmesh[i_rmin-1], rmesh[i_rmin])
    rmax = spopt.brentq(f_thrsh, rmesh[i_rmax], rmesh[i_rmax+1])

    return rmax - rmin


if __name__ == '__main__':
    main()
