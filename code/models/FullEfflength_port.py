"""
Port Sean's model code to Python, aiming for:
1. clarity
2. better speed and resolution

Main method is fefflen (contraction of FullEff_length), which takes
B0, eta2, mu, SNR parameters and outputs FWHMs in specified energy bands.

Aaron Tran
2014 August 6-8
"""

from __future__ import division

from datetime import datetime
import numpy as np
import scipy as sp
import scipy.optimize as spopt
import scipy.integrate as spint

# Recompile Fortran full model code and load Pacholczyk tables
#import fullmodel_recompile
#fullmodel_recompile.main()
import fullmodel
import snr_catalog as snrcat

# Constants
C_M = snrcat.SYNCH_CM
C_1 = snrcat.SYNCH_C1
A = snrcat.SYNCH_B  # Variable named A to avoid conflict w/ magnetic field B
NUKEV = snrcat.NUKEV

def main():
    """Currently used as a debugging testbed"""
    #B0 = raw_input('Enter B0 (Gauss): ')
    #eta2 = raw_input('Enter eta2: ')
    #mu = raw_input('Enter mu: ')
    B0 = 200e-6
    eta2 = 1.0
    mu = 1.0

    irmax, iradmax, ixmax = 200, 100, 500
    print 'irmax, iradmax, ixmax = {},{},{}'.format(irmax, iradmax, ixmax)

    kevs = np.array([0.7, 1., 2., 3., 4.5])
    rminarc = 20. * np.ones(len(kevs))

    # Tycho
    fwhms = fefflen(kevs, B0, eta2, mu, 4.52e8, 4.52e8/4, 1.077e19, 240.,
                    2*0.58+1, rminarc, True, irmax, iradmax, ixmax,
                    False, 0.001, 5e-6, 'fglists_mod.dat',
                    200, 50, 2000, False)
    # SN 1006 with damping
    #fwhms = fefflen(kevs, B0, eta2, mu, 5e8, 5e8/4., 2.96e19, 900.,
    #                2*0.6+1, rminarc, True, irmax, iradmax, ixmax,
    #                False, 0.001, 5e-6, 'fglists_mod.dat',
    #                200, 50, 2000, False)
    #print 'rminarc = {}'.format(rminarc)
    #for en, fwhm in zip(kevs, fwhms):
    #    print '{:0.2f} keV: {:0.5f}'.format(en, fwhm)

    return fwhms

# ==============================
# Main function to compute FWHMs
# ==============================

#@profile
def fefflen(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s, rminarc, icut, irmax,
            iradmax, ixmax, idamp, ab, Bmin, fgfname, itmax, inmax, irhomax,
            get_prfs):
    """Use numpy arrays for everything (kevs, rminarc)

    Note that irmax is not very important anymore.

    Input:
        kevs (np.ndarray)
        B0, eta2, mu (float), B0 in Gauss
        vs, v0, rs, rsarc, s (float)
        rminarc (np.ndarray)
        icut (bool, or 0/1)
        irmax, iradmax, ixmax (int)
            irmax = resolution on intensity profile (400, SN 1006)
            iradmax = resolution on tabulated e- distribution (100, SN 1006)
            ixmax = resolution for line of sight integration (500, SN 1006)
        idamp (bool), enable/disable damping
        ab, Bmin (float) damping scale length btwn [0,1], Bmin in Gauss
        fgfname (str) filename for 1-particle synchrotron emissivity table
        itmax, inmax, irhomax (int)
            itmax = resolution on finite, transformed Green's funct integral
                    for all cases (1000 default, SN 1006)
            inmax = resolution on infinite (then transformed) integral
                    for mu >=1 (100 default, SN 1006)
            irhomax = resolution in emistab.  This should exceed iradmax!
                      compute emissivity at irhomax pts, which is then summed
                      on line of sight (interpolating) to get intensity
                      irhomax < irmax should be okay.
        get_prfs (bool)
            if True, return computed intensity profile as a list of np.ndarray
    Output:
        np.ndarray of FWHMs for each entry in kevs
    """

    # Various derived constants
    rmin = 1 - rminarc / rsarc  # array of scaled r, range [0,1]
    # eta = \eta*(E_h)^(1-\mu), per section 4.2 of paper
    eta =  eta2 * (2.*NUKEV/(C_M*B0))**(-(mu-1)/2)
    if icut:
        Ecut = ( (8.3*(B0/(100e-6))**(-0.5)*(vs/1e8)*1.6021773)**(2/(1+mu))
                *(1./eta)**(1./(mu+1)) )  # I believe this is right.
        # Checked 2014 Sep 24 by Aaron
        # (1.602 erg)**(2/(1+mu)) balances (1/eta)**(1/(1+mu))
        # to give erg**1.  Everything looks okay.
    else:
        Ecut = 1e40

    # 1-particle synchrotron emissivity
    xex, fex = read_fglists(fgfname)
    fullmodel.readfglists(fgfname)  # !! important

    widths = []
    if get_prfs:
        rgrids = []
        prfs = []

    for en_nu, rmin_nu in zip(kevs, rmin):
        # Convert keV to s^{-1}
        nu = en_nu*NUKEV

        # Tabulate e- distr. for fixed nu over e- energy, radial coord
        # Be careful, Fortran code hardcoded for iradmax < 1000
        # and Pacholczyk table < 100 entries
        idamp_flag = 1 if idamp else 0  # Convert bool to int for Fortran
        radtab = np.linspace(rmin_nu, 1.0, num=iradmax)
        disttab = fullmodel.distr(np.empty(1000), iradmax, nu, rmin_nu,
                                  B0, Ecut, eta, mu, rs, v0, s, C_1,
                                  idamp_flag, ab, Bmin, itmax, inmax)
        disttab = disttab[:len(xex), :iradmax]  # Slice off zeros

        # Tabulate emissivity over radial coord
        rhotab = np.linspace(rmin_nu, 1., num=irhomax)
        emistab = emisx(rhotab, nu, B0, radtab, disttab, xex, fex,
                        idamp, ab, Bmin)

        # Make grids for intensity computation (x = line-of-sight coord)
        rmesh = np.linspace(rmin_nu, 1.0, num=irmax, endpoint=False)
        xmaxmesh = np.sqrt(1.-rmesh**2)  # Max x on line of sight at each r
        # Grid along line of sight for each r; shape is (irmax, ixmax)
        # [[0., ..., xmax_0], [0., ..., xmax_1], ... [0., ..., xmax_{irmax-1}]]
        xintv = np.linspace(0., 1., num=ixmax)  # Uniform x-coord sampling
        xgrid = xintv * xmaxmesh.reshape(-1,1)  # stackoverflow.com/a/16887425
        rgrid = np.tile(rmesh, (len(xintv), 1)).T  # r = constant on sightlines

        # Compute rho, emissivity at each point in grid (r,x)
        rhogrid = np.sqrt(xgrid**2 + rgrid**2)
        emisgrid = np.interp(rhogrid, rhotab, emistab)
        # then, integrate emis. along lines of sight
        intensity = spint.simps(emisgrid, xgrid)

        def get_intensity(r):
            """Non-vectorized intensity computation, for FWHM finding."""
            xmesh = np.linspace(0., np.sqrt(1.-r**2), num=ixmax)
            rhomesh = np.sqrt(xmesh**2 + r**2)
            emismesh = np.interp(rhomesh, rhotab, emistab)
            if r == 1:  # spint.simps fails for dx = 0.
                return 0.
            return spint.simps(emismesh, xmesh)

        # Finally, compute FWHMs precisely (use intensity grid as guide)
        w = fwhm(rmesh, intensity, get_intensity)
        widths.append(w)
        if get_prfs:
            rgrids.append(rmesh * rsarc)
            prfs.append(intensity)

    fwhms = np.array(widths) * rsarc

    if get_prfs:
        return fwhms, prfs, rgrids

    return fwhms


def read_fglists(fname):
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

def emisx(r, nu, B0, radtab, disttab, xex, fex,
          idamp, ab, Bmin):
    """Numerically integrate 1-particle synchrotron emissivity (emis1) over
    electron distribution to obtain full emissivity j_\nu at r (array)

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
    # Would be nice to use np.interp, but not sure if possible
    ind = np.searchsorted(radtab, r)
    slp = (disttab[:,ind] - disttab[:,ind-1]) / (radtab[ind] - radtab[ind-1])
    dist = disttab[:,ind-1] + slp * (r - radtab[ind-1])
    # dist = slice of disttab at fixed x = 1-r, linearly interpolated

    xex = np.tile(xex, (len(r),1)).T  # e- energy
    fex = np.tile(fex, (len(r),1)).T  # e- single-particle emissivity

    integrand = fex * dist / (xex**1.5)

    # Damping: scale emissivity by sqrt(nu*B(r))
    if idamp:
        B = Bmin + (B0 - Bmin)*np.exp(-1*(1.0 - r) / ab)
    else:
        B = B0

    # Skip constant prefactors
    return spint.trapz(integrand, xex, axis=0) * np.sqrt(nu*B)


# ===============================
# Electron distribution functions
# ===============================

# NOTE distr(...) not currently in use!  Calling Fortran code instead

#def distr(Ec, rc, B, Ecut, eta, mu, rs, v, s):
#    """Electron distribution, from Rettig & Pohl and Lerche & Schlickeiser
#
#    ASSUMED: E, r are 1-D arrays
#    spint.trapz integrates along last axis by default, conveniently
#
#    Input:
#        E = particle energy, r = radial position (scaled to [0,1])
#        B = magnetic field
#        Ecut = cut-off energy of injected electrons
#        eta = \eta*(E_h)^(1-\mu), from the paper
#        mu = diffusion coefficient exponent
#        rs, v, s = shock radius, plasma velocity, injected e- spectral index
#    """
#
#    D0 = eta*2.083e19/B  # eta (E_h)^(1-mu) * C_d / B0
#    a = 1.57e-3  # This is b in the paper
#                 # Synchrotron constant ... TODO declare at module level?
#    k0 = 1.  # Normalization, Q_0
#    b0 = a * B**2  # b*B^2 in paper
#    alpha = 1. # TODO Not sure what this is, scaling factor?
#
#    def get_nums(E, r):
#        """Compute derived parameters over 3-D arrays from np.meshgrid"""
#        tau = 1. / (b0 * E)  # Synchrotron cooling time
#        lad = v * tau  # Advective lengthscale
#        ldiff = np.sqrt(D0*E**mu*tau)  # Diffusive lengthscale
#        x = (rs - r*rs)/alpha  # Convert scaled r to x in cm
#        return tau, lad, ldiff, x
#
#    # Functions for each of 4 cases
#
#    # TODO this one is not yet functional
#    def distr_adv(Ec, rc, B, Ecut, eta, mu, rs, v, s):
#        """Case: l_ad / l_diff >> 1"""
#        efinv = a/v * B**2 * (rs - r*rs)
#        ef = 1. / efinv
#        en0 = E / (1.0 - E/ef)
#        argexp = en0/Ecut
#        # argexp = E / Ecut / (1 - a*B**2*E * (rs-r*rs) / v)
#        Xi = 1 - E/ef  # Only used for if-else structure below
#        if Xi > 0:
#            return 1./en0**s * (en0/E)**2 / np.exp(argexp) * 8e-9
#        else:
#            return 0.
#
#
#    def distr_mlt1(Ec, rc, B, Ecut, eta, mu, rs, v, s):
#        """Case: mu < 1
#        Integrate over t, where n = 1 - t^2
#        """
#        # Build 3-D arrays
#        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
#        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
#        n = 1 - t**2
#        tau, lad, ldiff, x = get_nums(E, r)
#
#        argexp = (n**(-1/(1-mu)) * E/Ecut +
#                      (lad/alpha * (1 - n**(1/(1-mu)) ) - x)**2 /
#                      (4*ldiff**2*(1-n)/alpha) * (1-mu) )
#        argexp[n==1.] = 1e35  # Prevent blowup
#        integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)
#
#        integral = spint.trapz(integrand, t)  # TODO can do better?
#        E_s = E[:,:,0]  # Slice as integral is a 2-D array
#        return integral * ( k0 /2 /np.sqrt(np.pi*alpha*b0*D0*(1-mu)) *
#                            E_s**(-1*(mu/2+1/2+s)) )
#
#
#    def distr_rpohl(Ec, rc, B, Ecut, eta, mu, rs, v, s):
#        """Case: mu = 1
#        Integrate over t in [0,1] (n=e^(t^2) in [1,e]), then
#        integrate over y in [0,1] (n=mess in [e, infty)).
#        """
#        # First integral, build 3-D arrays
#        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
#        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
#        n = np.exp(t**2)
#        tau, lad, ldiff, x = get_nums(E, r)
#
#        argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
#                              (4*ldiff**2/alpha*np.log(n)) )
#        argexp[n==1.] = 1e35  # Prevent blowup
#        integrand = 2*np.exp((1-s)*t**2) / np.exp(argexp)
#        integral_left = spint.trapz(integrand, t)  # TODO can do better?
#
#        # Second integral, build 3-D arrays
#        yc = np.linspace(0., 1., num=100, endpoint=False)  # TODO arbitrary
#        E, r, y = np.meshgrid(Ec, rc, yc, indexing='ij')
#        q = 2.  # q, nmin just for change of variables
#        nmin = np.exp(1.)
#        n = y / (1 - y**2)**q + nmin
#        tau, lad, ldiff, x = get_nums(E, r)
#
#        argexp = n*E/Ecut + ( (lad/alpha*(1-1/n) - x)**2/
#                              (4*ldiff**2/alpha*np.log(n)) )  # same
#        integrand = ( n**(-1*s) / np.sqrt(np.log(n)) / np.exp(argexp) *
#                      ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
#        integral_right = spint.trapz(integrand, y)  # TODO can do better?
#
#        # Final integral
#        integral = integral_left + integral_right
#        E_s = E[:,:,0]  # Slice as integral is 2-D array
#        return integral * k0 /2 /np.sqrt(np.pi*alpha*b0*D0) * E_s**(-1*(1+s))
#
#
#    def distr_mgt1(Ec, rc, B, Ecut, eta, mu, rs, v, s):
#        """Case: mu > 1
#        Integrate over t in [0, 1] (n = 1+t^2 in [1, 2]), then
#        integrate over y in [0, 1] (n = mess in [2, infty)).
#        """
#        # First integral, build 3-D arrays
#        tc = np.linspace(0., 1., num=1000)  # TODO arbitrary
#        E, r, t = np.meshgrid(Ec, rc, tc, indexing='ij')
#        n = 1 + t**2
#        tau, lad, ldiff, x = get_nums(E, r)
#
#        argexp = n**(-1/(1-mu)) * E/Ecut + ( (lad/alpha*(1-n**(1/(1-mu)))-x)**2 /
#                                             (4*ldiff**2*(1-n)/alpha) * (1-mu) )
#        argexp[n==1.] = 1e35
#        integrand = 2*n**((s+mu-2)/(1-mu)) / np.exp(argexp)
#        integral_left = spint.trapz(integrand, t)  # TODO can do better?
#
#        # Second integral, build 3-D arrays
#        yc = np.linspace(0., 1., num=100, endpoint=False)  # TODO arbitrary
#        E, r, y = np.meshgrid(Ec, rc, yc, indexing='ij')
#        q = 2.
#        nmin = 2.
#        n = y / (1 - y**2)**q + nmin
#        tau, lad, ldiff, x = get_nums(E, r)
#
#        argexp = ( n**(-1/(1-mu)) * E/Ecut +
#                   (lad/alpha*(1-n**(1/(1-mu))) - x)**2 /
#                   (4*ldiff**2*(1-n)/alpha) * (1-mu) )
#        integrand = ( n**((s+mu-2)/(1-mu)) / np.sqrt(n-1) / np.exp(argexp) *
#                      ((2*q-1)*y**2+1) / (1-y**2)**(q+1) )
#        integral_right = spint.trapz(integrand, y)  # TODO can do better?
#
#        integral = integral_left + integral_right
#        E_s = E[:,:,0]  # Slice as integral is a 2-D array
#        return integral * ( k0 /2 /np.sqrt(np.pi*alpha*b0*D0*(mu-1)) *
#                            E_s**(-1*(mu/2+1/2+s)) )
#
#    # Need to address edge case by vectorizing... or something.
#    # Construct 2-D grid of E, r
#
#    # lad/ldiff ~ sqrt(v*lad) / dsqrt(D) ~ peclet number!
#    #if lad/ldiff > 30:
#    #    return distr_adv(E, r, B, Ecut, eta, mu, rs, v, s)
#
#    if mu > 1:
#        return distr_mgt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
#    elif mu < 1:
#        return distr_mlt1(Ec, rc, B, Ecut, eta, mu, rs, v, s)
#    else:
#        return distr_rpohl(Ec, rc, B, Ecut, eta, mu, rs, v, s)


# ============================
# Find FWHM in intensity graph
# ============================

def fwhm(rmesh, intensity, f_int):
    """Compute fwhm precisely from starting grid.
    For extreme values of B0/eta2 (extremely smeared, no peak exists),
    code prints error messages and returns FWHM = 1 (max possible value)"""
    # Tolerance for finding intensity max position
    eps2 = np.finfo(float).eps * 2  # Tolerance, 1 ~= 1+eps2/2.

    # Compute half max from grid initial guess by bracketing intensity max
    idxmax = np.argmax(intensity)
    if idxmax == 0:
        print 'ERROR: intensity max not found (rminarc too small)'
        return 1  # Can't search r < rmin, no disttab/emisttab values computed
    elif idxmax == len(intensity)-1:
        #print 'Warning: intensity max not resolved in grid on right; searching'
        rpk_a = rmesh[-2]  # CANNOT be rmesh[-1], to find max
        rpk_b = search_crossing(rmesh[-1], 1., lambda r:f_int(r)-intensity[-1],
                                eps2)  # Find r s.t. f_int(r) < intensity[-1]
        if rpk_b is None:  # This should never happen, honestly
            print 'ERROR: intensity max not found (stalled on right) (?!)'
            return 1
    else:
        rpk_a = rmesh[idxmax - 1]
        rpk_b = rmesh[idxmax + 1]

    # option 'xatol' requires SciPy 0.14.0
    # 'xatol' -- absolute error in res.x acceptable for convergence
    # (as res.x is order 1, eps2 should be appropriate)
    res = spopt.minimize_scalar(lambda x: -1*f_int(x), method='bounded',
                                bounds=(rpk_a, rpk_b), options={'xatol':eps2})
    rpk = res.x
    pk = f_int(rpk)
    halfpk = 0.5 * pk
    def f_thrsh(r):  # For rootfinding
        return f_int(r) - halfpk

    # Right (upstream) FWHM crossing -- do not search on grid
    # (grid is too coarse upstream of rim, will not find FWHM;
    #  spurious crossings occur if profile not monotonic, e.g. w/ damping)
    rmax = spopt.bisect(f_thrsh, rpk, 1.)

    # Left (downstream) FWHM crossing -- find bracketing r (position) values
    # requires that rminarc is large enough so that grid contains crossing
    cross = np.diff(np.sign(intensity - halfpk))
    inds_rmin = np.where(cross > 0)[0]  # Left (neg to pos)

    if inds_rmin.size == 0:  # No left crossing found
        print ('ERROR: FWHM edge (rmin) not found '
               '(rminarc too small or peak FWHM cannot be resolved in profile)')
        return 1.  # Can't search r < rmin, no disttab/emisttab values computed
    else:
        rmin_a = rmesh[inds_rmin[-1]]  # Crossing closest to peak (largest r)
        rmin_b = rmesh[inds_rmin[-1] + 1]

    rmin = spopt.bisect(f_thrsh, rmin_a, rmin_b)

    return rmax - rmin


def search_crossing(r_init, r_lim, f, eps):
    """Find r s.t. f(r) < 0, up to r_limit by binary search
    (i.e., assumes monotonicity), with f(r_init) > 0.
    If no crossing is found, return None.
    eps should be positive.
    """
    r = r_init + (r_lim - r_init) / 2.  # Increment r right away to avoid
    while f(r) >= 0:                    # floating point error in f(r_init)
        r += (r_lim - r) / 2.
        if abs(r_lim - r) < eps:  # Breaking criterion
            return None
    return r


if __name__ == '__main__':
    main()
