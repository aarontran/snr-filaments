"""
Catalog of SNR physical parameters and synchrotron constants

Aaron Tran, Summer 2014
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
                    # //cv.nrao.edu/course/astr534/SynchrotronSpectrum.html
                    # anyways this is all one giant TODO that should be typed
                    # up and thrown into the code deep review PDF

CD = 2.083e19  # Bohm diffusion const; Cd = c/(3e)
BETA = 4.6  # Projection factor for exponential emissivity (Ballet, 2006)

NUKEV = 2.417989e17  # Frequency of 1 keV photon

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
    """Physical parameters, model/fit settings for Tycho's SNR
    WARNING: if darc is changed, must change vs as well!
    """
    k = SupernovaRemnant('Kepler')

    kepler.dkpc = 3.3  # Distance to remnant in kpc (!) TODO verify
                       # 4kpc from Vink (2008) kinematic study
                       # Could be up to 7 kpc?!?!?! (Patnaude et al., 2012)
                       # 3.3kpc from Katsuda et al. (2008)
    kepler.rsarc = 180 # Shock radius (arcsec) from Green's SNR catalog
    kepler.s = 2.28  # e- spectrum index, 2.3 = 2*0.65 + 1;
                     # 0.64 = radio spectral index, Green's SNR catalog
                     # TODO get a source on this?!

    kepler.vs = 4.2e8  # Shock velocity, cm/s
                       # TODO verify, 4.2e8 from Vink (2008) kinematic study
                       # Katsuda et al assumed 1.6e8 (?!)
                       # if shock velocity is truly higher, then distance
                       # estimate of Katsuda et al. should be revised up...
    kepler.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    kepler.icut = True
    #kepler.rminarc = 20  # Default rminarc, arcsec
    kepler.irmax = 400  # Currently unchanged from SN 1006
    kepler.iradmax = 100
    kepler.ixmax = 500

    # Fitting default initial guesses, bounds
    #kepler.par_init = {'mu': 1.0, 'eta2': 1.0, 'B0': 300e-6}
    #kepler.par_lims = {'mu': (0., 2.),
    #                  'eta2': (1e-16, 1e5),  # model code div by zero on eta2=0
    #                  'B0': (1e-6, 1e-2)}

    return kepler

def make_tycho():
    """Physical parameters, model/fit settings for Tycho's SNR
    WARNING: if darc is changed, must change vs as well!
    """
    tycho = SupernovaRemnant('Tycho')

    tycho.dkpc = 3.0  # Distance to remnant in kpc (!)
    tycho.rsarc = 240  # Shock radius (arcsec) from Green's SNR catalog
    tycho.s = 2.3  # e- spectrum index, 2.3 = 2*0.65 + 1;
                   # 0.65 = radio spectral index
                   # TODO is this right?  One src. is Eriksen et al. 2011
                   # which cites Kothes et al. 2006, for alpha=0.65
                   # for the ENTIRE remnant (mean).
                   # Katz-Stone et al. (2000) suggest 0.5 for various
                   # "filaments" within the remnant?
                   # And, Green's SNR catalog gives 0.58, not 0.65.  DAMMIT.

    tycho.vs = 3.6e8 * tycho.dkpc/2.3  # Average-ish shock velocity, cm/s
    # Mean of velocs from nonthermal regions, Williams et al. 2013
    tycho.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    tycho.icut = True
    tycho.rminarc = 20  # Default rminarc, arcsec
    tycho.irmax = 400  # Currently unchanged from SN 1006
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

    sn1006.vs = 5e8  # Shock velocity, cm/s (Satoru et al., 2009, 2013)
    sn1006.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    # Model grid settings, resolutions
    sn1006.icut = True
    sn1006.rminarc = 60  # Default rminarc, arcsec
    sn1006.irmax = 400
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
