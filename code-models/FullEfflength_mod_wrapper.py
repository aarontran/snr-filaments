"""
Code to f2py-ize FullEfflength_mod.f
and attempt to run fits...

FullEfflength_mod.f modified appropriately
(just revert to older version as needed)

Based on code-profiles/fwhm-process.ipynb setup for the simpler code
(just a prototype).  I will probably use this code for both models, to
consolidate global constants + fit work all in one place.

Aaron Tran
2014 July 21
"""

import lmfit
import matplotlib.pyplot as plt
import numpy as np
from numpy import f2py
import scipy as sp
from scipy import optimize

with open('FullEfflength_mod.f', 'r') as f:
    fsource = ''.join(list(f))
np.f2py.compile(fsource, modulename='fullmodel', verbose=1)
import fullmodel as fm
fm.readfglists()

def main():
    """Docstring me"""


# =========================
# Run fits for whatever SNR
# =========================

def apply_fit():
    """Just testing..."""
    kevs = np.array([0.7, 1., 2.])
    #data4 = np.array([29, 23.9, 16.6])
    #eps4 = np.array([.9, .39, .45])
    data5 = np.array([33.75, 27.2, 24.75 ])
    eps_data5 = np.array([2.37,.62,.61])

    data, eps = data5, eps_data5

    #p = lmfit.Parameters()
    #p.add('B0', value=1e-4, min=1e-6, max=1e-2)
    #p.add('eta2', value=1, vary=False) #min=0, max=1e4)
    #p.add('mu', value=1, vary=False)

    #res = lmfit.minimize(width_cont, p,
    #               args=(kevs, data, eps),
    #               method='leastsq')

    #print lmfit.printfuncs.report_fit(res.params)

    f = lambda b: width_cont(b, kevs, data, eps)

    popt, pcov = sp.optimize.leastsq(f, np.array([1e-4]), epsfcn=1e-6)

    #en_plt = np.linspace(0.6, 2.1, 10)
    #plt.errorbar(kevs, data, eps, fmt='bo')
    #plt.plot(en_plt, width_cont(p, en_plt), '-k')
    #plt.xlabel('Observation energy (keV)')
    #plt.ylabel('FWHM (arcsec)')
    #plt.show()

    #return res
    return popt, pcov


def params_sn1006():
    """Configure parameters for SN1006
    output: 
    """
    compratio = 4.0  # Compression ratio from Rankine-Hugoniot

    p = lmfit.Parameters()

    p.add('B0', value=1e-4, min=1e-6, max=1e-2)
    p.add('eta2', value=1, min=0, max=1e4)

    p.add('mu', value=1, vary=False)

    #p.add('vs', value=5e8, vary=False)  # Shock velocity (cm/s)
    #p.add('v0', value=5e8/compratio, vary=False)  # Plasma velocity (cm/s)
    #p.add('rs', value=2.96e19, vary=False)  # Shock radius (cm)
    #p.add('rsarc', value=900, vary=False)  # rs in arcsec, Green's SNR catalog
    #p.add('s', value=2.2, vary=False)  # 2*0.6 + 1, 0.6 = radio spectral index

    #p.add('rminarc', value=60, vary=False)
    #p.add('icut', value=0, vary=False)
    #p.add('irmax', value=400, vary=False)
    #p.add('iradmax', value=100, vary=False)
    #p.add('ixmax', value=500, vary=False)

    return p


# =====================
# Functions for fitting
# =====================

#def width_cont(params, kevs, dat, eps):
def width_cont(B0, kevs, dat, eps):
    """Width function to fit; uses numerical output from full model
    Wrapper for fm.fullefflengthsub -- feeds in the correct input,
    outputs the desired output (that is, width values in arcsec)
    """
    # Call signature, for reference (KEEP THIS UPDATED)
    # Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu,
    #                  vs, v0, rs, rsarc, s, rminarc, icut,
    #                  irmax, iradmax, ixmax)

    # f2py sets widtharc as output (from cf2py directive) and inumax optional
    # Call signature is, then:
    # fullefflengthsub(kevs, b0, eta2, mu, vs, v0, rs, rsarc, s, rminarc
    #                  icut, irmax, iradmax, ixmax, [inumax])

    #B0 = params['B0'].value
    #eta2 = params['eta2'].value
    #mu = params['mu'].value
    eta2 = 1; mu = 1

    vs = 5e8
    v0 = 5e8/4.
    rs = 2.96e19
    rsarc = 900
    s = 2.2

    rminarc = 60
    icut = 0
    irmax = 400
    iradmax = 100
    ixmax = 500

    #vs = params['vs'].value
    #v0 = params['v0'].value
    #rs = params['rs'].value
    #rsarc = params['rsarc'].value
    #s = params['s'].value

    #rminarc = params['rminarc'].value
    #icut = params['icut'].value
    #irmax = params['irmax'].value
    #iradmax = params['iradmax'].value
    #ixmax = params['ixmax'].value

    print 'function call, b0 = {} microG'.format(B0*1e6)

    return (fm.fullefflengthsub(kevs, B0, eta2, mu, vs, v0, rs, rsarc, s,
                                rminarc, icut, irmax, iradmax, ixmax)
            - dat) / eps


def width_dump(params, kevs):
    """Width function to fit, from Parizot et al."""
    pass


def objectify(f):
    """Generate objective function from model function + data with weights
    Yes, it's named objectify.
    
    Inputs:
        f: model function for data, signature f(pars, x)
    Inputs to output function:
        data: measured data to fit
        eps_data: errors on measured data (used as weights)
        *args: necessary arguments to f; for a simple 1-D curve fit, this
        will include x-coordinates (abscissae) corresponding to values of data
    Output:
        objective function to minimize in least-squares sense, or otherwise,
        using lmfit.minimize or scipy.optimize.leastsq
    """
    return lambda pars, x, data, eps_data: (f(pars, x) - data)/eps_data


if __name__ == '__main__':
    main()
