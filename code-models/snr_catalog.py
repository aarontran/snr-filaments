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

def make_tycho():
    """Physical parameters for Tycho's SNR
    WARNING: if darc is changed, must change vs as well!
    """
    tycho = SupernovaRemnant('Tycho')
    
    tycho.dkpc = 3.0  # Distance to remnant in kpc (!)
    tycho.rsarc = 240  # Shock radius (arcsec) from Green's SNR catalog
    tycho.s = 2.3  # e- spectrum index, 2.3 = 2*0.65 + 1;
                   # 0.65 = radio spectral index

    tycho.vs = 3.6e8 * tycho.dkpc/2.3  # Shock velocity, cm/s
    # Mean of velocs from nonthermal regions, Williams et al. 2013
    tycho.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    tycho.rminarc = 12  # Default rminarc, arcsec
    #tycho.kevs = [0.7, 1., 2., 3., 4.5]  # Photon energies
    #tycho.inds = None  # None implies use all bands

    return tycho


def make_SN1006():
    """Physical parameters for SN 1006"""
    sn1006 = SupernovaRemnant('SN 1006')

    sn1006.dkpc = 2.2  # Distance to remnant in kpc (!!)
    sn1006.rsarc = 900  # Shock radius (arcsec) from Green's SNR catalog
    sn1006.s = 2.2  # e- spectrum index, 2.2 = 2*0.6 + 1;
                    # 0.6 = radio spectral index

    sn1006.vs = 5e8  # Shock velocity, cm/s (Satoru et al., 2009, 2013)
    sn1006.cratio = 4.0  # Compression ratio, strong adiabatic shock (R-H)

    sn1006.rminarc = 60  # Default rminarc, arcsec
    #sn1006.kevs = [0.7, 1., 2.]  # Photon energies
    #sn1006.inds = None  # None implies use all bands

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


def main():
    pass


if __name__ == '__main__':
    pass
