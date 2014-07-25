"""
Catalog of SNR physical parameters and synchrotron constants

Aaron Tran, Summer 2014
2014 July 22
"""

# TODO verify that constants are correct!
# Constants are also defined in Fortran code, check those as well

# All constants in CGS (Gaussian) units
SYNCH_B = 1.57e-3  # Synchrotron cooling time const
SYNCH_CM = 1.82e18  # Synchrotron characteristic freq const

CD = 2.083e19  # Bohm diffusion const; Cd = c/(3e)
BETA = 4.6  # Projection factor for exponential emissivity (Ballet, 2006)

NUKEV = 2.417989e17  # Frequency of 1 keV photon


class SupernovaRemnant(object):
    def __init__(self, name):
        self.name = name


def make_tycho():
    """Physical parameters for Tycho's SNR"""
    tycho = SupernovaRemnant('Tycho')
    
    tycho.vs = 3e8  # Approx. shock velocity cm/s (cf. Williams et al., 2013)
                    # Shock velocity is variable around SNR, account for that!
    tycho.v0 = tycho.vs/4.0  # Plasma velocity, cm/s (compratio = 4)
    tycho.rsarc = 240  # Shock radius (arcsec) from Green's SNR catalog
    tycho.rs = 1.077e19  # Shock radius (cm), tan(240 arcsec) * 3 kpc
    tycho.s = 2.3  # e- spectrum index, 2.3 = 2*0.65 + 1; 0.65 = radio spectral index

    return tycho


def make_SN1006():
    """Physical parameters for SN 1006"""
    sn1006 = SupernovaRemnant('SN 1006')

    sn1006.vs = 5e8  # Shock velocity, cm/s (Satoru et al., 2009, 2013)
    sn1006.v0 = sn1006.vs/4.0  # Plasma velocity, cm/s (compratio = 4)
    sn1006.rsarc = 900  # Shock radius (arcsec) from Green's SNR catalog
    sn1006.rs = 2.96e19  # Shock radius (cm), tan(900 arcsec) * 2.2 kpc
    sn1006.s = 2.2  # e- spectrum index, 2.2 = 2*0.6 + 1; 0.6 = radio spectral index

    return sn1006


def main():
    pass

if __name__ == '__main__':
    pass
