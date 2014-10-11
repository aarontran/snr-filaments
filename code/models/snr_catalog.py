"""
Catalog of SNR physical parameters and synchrotron constants

Aaron Tran
2014 July 22
"""

from __future__ import division
import numpy as np

# TODO verify that constants are correct!
# Constants are also defined in Fortran code, check those as well
# All constants in CGS (Gaussian) units

SYNCH_B = 1.57e-3  # Synchrotron cooling time const
SYNCH_CM = 1.82e18  # Synchrotron characteristic freq const
SYNCH_C1 = 6.27e18  # for synchrotron emissivity (eq'n 20; see Pacholczyk)
                    # I believe C_1 = 3e/(4\pi m^3 c^5), in CGS units
                    # e = e- charge, m = e- mass, c = 3e10 cm/s
                    # i.e., prefactor for critical freq. as in
                    # cv.nrao.edu/course/astr534/SynchrotronSpectrum.html
                    # anyways this is all one giant TODO that should be typed
                    # up and thrown into the code deep review PDF

CD = 2.083e19  # Bohm diffusion const; Cd = c/(3e) (confirmed)
BETA = 4.6  # Projection factor for exponential emissivity (Ballet, 2006)

NUKEV = 2.417989e17  # Frequency of 1 keV photon (confirmed)

# ================================
# Some convenient unit conversions
# ================================

def arcsec2rad(x):
    """Convert arcsec to radians"""
    return x / 3600. * np.pi/180.

def kpc2cm(x):
    """Convert kiloparsec to cm.
    1 parsec = 3.085678e16 m (NIST guide to SI units, 2008)"""
    return x * 1e3 * 3.085678e16 * 1e2


# ================================
# Set SNR physical parameters here
# ================================

def make_kepler():
    """Physical parameters, model/fit settings for Kepler's SNR
    WARNING: if dkpc is changed, must change vs as well!
    """
    kepler = SupernovaRemnant('Kepler')

    kepler.dkpc = 5    # Distance to remnant in kpc
                       # Many conflicting results -- X-ray expansion gives
                       # numbers of ~3.3 kpc (Katsuda et al., 2008)
                       # Vink (2008) uses a funny maximum likelihood method
                       # and gets 4 kpc, but favors much larger distances
                       # Chiotellis et al. (2012) and Patnaude et al. (2012)
                       # favor >6,7 kpc for standard (1 M_56Ni or w/e)
                       # energetics, but ~4-6 kpc consistent w/ subenergetic
                       # explosion.
                       # TeV non-detection by HESS favors > 6.4 kpc
                       # (Aharonian et al., 2008)
                       # Finally, radio HI absoprtion gives [4.8, 6.4] kpc
                       # (Reynoso and Goss, 1999)
                       # And Balmer filament measurements give ~4 kpc or so...
                       # (Sankrit et al., 2005)
                       # Choice of 5 kpc also taken by Toledo-Roy et al. (2014)

    kepler.rsarc = 114 # Shock radius (arcsec) = 1.9 arcmin, compromise between
                       # ~2 arcmin at ears and 1.8 arcmin at rest of remnant
                       # Green's estimate of 1.5 arcmin is too small
                       # Chiotellis et al. (2012) use ~1.7-1.8 arcmin

    kepler.s = 2.28  # e- spectrum index, 2.28 = 2*0.64 + 1;
                     # 0.64 = radio spectral index, Green's SNR catalog
                     # consistent w/ DeLaney et al. (2002)

    kepler.vs = 4.71e8 * kepler.dkpc/5 # Shock velocity, cm/s (with dkpc = 5)
                        # mean measurements of Reg-4/5, Katsuda et al. (2008)
                        # scaled to 5 kpc.  Vink (2008) give 5.25e8
                        # but within error the values are barely consistent

    kepler.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    kepler.icut = True
    kepler.rminarc = 20  # Default rminarc, arcsec (remember w/ adaptive
                         # rminarc calculation, this isn't too important)
    kepler.irmax = 200  # Currently unchanged from SN 1006
    kepler.iradmax = 100
    kepler.ixmax = 500

    # Fitting default initial guesses, bounds
    kepler.par_init = {'mu': 1.0, 'eta2': 1.0, 'B0': 250e-6}
    kepler.par_lims = {'mu': (0., 2.),
                       'eta2': (1e-16, 1e5),  # model code div by zero on eta2=0
                       'B0': (1e-6, 1e-2)}

    return kepler

def make_tycho():
    """Physical parameters, model/fit settings for Tycho's SNR
    WARNING: if dkpc is changed, must change vs as well!
    """
    tycho = SupernovaRemnant('Tycho')

    tycho.dkpc = 3.0  # Distance to remnant in kpc (!)
    tycho.rsarc = 240  # Shock radius (arcsec) from Green's SNR catalog
                       # This is about right from image inspection.
    tycho.s = 2.16 # e- spectrum index, 2.16 = 2*0.58 + 1
                   # 0.58 = radio spectral index from Sun et al. (2011)
                   # also listed in Green's SNR catalog.
                   # Before 2014 Oct 10 We used 0.65 from Kothes et al. (2006)

    tycho.vs = 3.6e8 * tycho.dkpc/2.3  # Average-ish shock velocity, cm/s
    # Mean of velocs from nonthermal regions, Williams et al. 2013
    tycho.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    tycho.icut = True
    tycho.rminarc = 20  # Default rminarc, arcsec
    tycho.irmax = 200
    tycho.iradmax = 100
    tycho.ixmax = 500

    # Fitting default initial guesses, bounds
    tycho.par_init = {'mu': 1.0, 'eta2': 1.0, 'B0': 300e-6}
    tycho.par_lims = {'mu': (0., 2.),
                      'eta2': (1e-16, 1e5),  # model code div by zero on eta2=0
                      'B0': (1e-6, 1e-2)}

    return tycho


def make_SN1006():
    """Physical parameters, model/fit settings for SN 1006"""
    sn1006 = SupernovaRemnant('SN1006')

    sn1006.dkpc = 2.2  # Distance to remnant in kpc (!!)
    sn1006.rsarc = 900  # Shock radius (arcsec) from Green's SNR catalog
    sn1006.s = 2.2  # e- spectrum index, 2.2 = 2*0.6 + 1;
                    # 0.6 = radio spectral index

    sn1006.vs = 5e8  * sn1006.dkpc / 2.2  # Shock velocity, cm/s
                                          # (Satoru et al., 2009, 2013)
    sn1006.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    sn1006.icut = True
    sn1006.rminarc = 60  # Default rminarc, arcsec
    sn1006.irmax = 200
    sn1006.iradmax = 100
    sn1006.ixmax = 500

    # Fitting default initial guesses, bounds
    sn1006.par_init = {'mu': 1.0, 'eta2': 1.0, 'B0': 150e-6}
    sn1006.par_lims = {'mu': (0., 2.),
                       'eta2': (1e-16, 1e5),  # model divides by zero on eta2=0
                       'B0': (1e-6, 1e-2)}

    return sn1006


# =======================================
# Define functions for derived parameters
# =======================================

class SupernovaRemnant(object):
    def __init__(self, name):
        self.name = name

    def v0(self):
        """Plasma velocity (cm/s)"""
        return self.vs / self.cratio

    def rs(self):
        """Shock radius (cm)"""
        return arcsec2rad(self.rsarc) * kpc2cm(self.dkpc)

    def config_log(self):
        """Dump information about current settings in printable string"""
        log = ['Properties:']
        for v in vars(self).items():
            log.append('{} = {}'.format(*v))
        log.append('Derived parameters:')
        log.append('v0 = {}'.format(self.v0()))
        log.append('rs = {}'.format(self.rs()))
        return '\n'.join(log)


if __name__ == '__main__':
    pass
