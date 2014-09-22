"""
Initialize and populate regdict objects for profile fitting, data storage
Designed to feed plot output into iPython notebooks for interactive work.

regdict is a dict keyed by energy band labels, w/ band dictionaries as values
band dictionaries (bdicts) store data, FWHMs, other useful information

rdict keys:
    '0.7-1kev', '1-2kev', '2-7kev', etc... tied to bdicts
bdict keys (when populated):
    'data' (x, y, y_err, c)
    'wtype'
    'wlen'
    'fwhm', 'lims', 'errs'
    'f', 'pars'

The general approach is to glob the region dictionaries, then iterate over the
regions and energy bands.  The profile data, wtype/wlen, fwhms, and f/pars
should provide enough for subsequent analysis and plotting.

Aaron Tran
September 2014
"""

import cPickle as pickle
from datetime import datetime
from glob import glob
import inspect
import json
import matplotlib.pyplot as plt
import numpy as np
import re

import ffwhm
import fsmooth
from json_utils import NumpyJSONEncoder
import regdict_utils as rdu
from regparse import check_dir

from fplot import fplot

# ==================================
# Methods to run in IPython notebook
# ==================================

def main_check(inroot, labels, wlen, wtype):
    """Main method to review fit smoothing before running fits
    Unlike other scripts, should be run in an IPython notebook
    """
    for rdict, n_reg in regdict_init(inroot, labels, wlen, wtype):
        print 'Region: {}'.format(n_reg)
        plt.figure(figsize=(5*len(labels), 4))
        for lab, i in zip(labels, xrange(len(labels))):
            plt.subplot(1, len(labels), i+1)
            rdu.plot_smoothing_cut(rdict[lab], plt.gca())
            plt.title(lab)
        plt.show()

def main_fwhm(inroot, outroot, labels, fitter, blacklist, wlen, wtype,
              save=True, **kwargs):
    """Main method to compute FWHMs and generate lots of output
    Unlike other scripts, should be run in an IPython notebook.

    Available **kwargs are:
        dx (float, arcsec)
        xconst (float, arcsec)
        cap (bool)
        sub_bkg (bool)
    """
    spec_cuts = []

    for rdict, n_reg in regdict_init(inroot, labels, wlen, wtype):
        blist_reg = blacklist[n_reg] if n_reg in blacklist else None
        rdict = get_fwhm_fits(rdict, n_reg, labels, fitter, blist_reg,**kwargs)

        if save:
            x_min, x_btw, x_max = get_cuts(rdict)
            spec_cuts.append([n_reg, x_min, x_btw, x_max])
            regdict_dump(rdict, n_reg, outroot)

    if save:
        speccut_dump(spec_cuts, outroot)

# ========================================
# Methods to create, populate region dicts
# ========================================

def regdict_init(inroot, labels, wlen, wtype):
    """Generate region dictionary for each region in profile data
    Depends on many settings/filenames in prep_profile_fit.py source!

    Input:
        inroot: filestem for processed profiles from prep_profile_fit.py
        labels: list of energy-band labels, must be in order!
                (at least, order must match input to prep_profile_fit.py)
        wlen, wtype: smoothing window params used to compute fit cuts
    Yield:
        a 2-tuple of
        1. region dict (keys: labels) of band dicts w/ 'data', 'wtype', 'wlen'
        2. region number
    """
    nregs = len(glob('{}_*_band_{}.npz'.format(inroot, labels[0])))
    for n in xrange(nregs):
        cuts = np.load('{}_fit_cuts.npz'.format(inroot), 'r')['cuts']
        cuts = cuts[cuts[:,0] == n+1][0][1:]  # Get cuts for region number n

        rdict = {}
        for lab, c in zip(labels, cuts):
            rdict[lab] = {}
            prfs = np.load('{}_{:02d}_band_{}.npz'.format(inroot, n+1, lab))
            rdict[lab]['data'] = prfs['x'], prfs['y'], prfs['y_err'], c
            rdict[lab]['wtype'] = wtype
            rdict[lab]['wlen'] = wlen

        yield rdict, n+1

def get_fwhm_fits(rdict, n_reg, labels, fitter, blist_reg=None, **kwargs):
    """Compute and store FWHM fits for all rdict bands
    Intended for interactive IPython notebook usage
    Input:
        rdict, n_reg: initialized region dict, region number
        fitter: function passed to bdict_fwhm_fit for profile fitting
        blist_reg: list of blacklisted (bad) energy bands, else None
        Other **kwargs (cap, sub_bkg, dx, xconst) sent to bdict_fwhm_fit
    Output:
        region dict w/data populated for each energy band
    """
    fig, axes = plt.subplots(1, len(rdict), figsize=(5*len(rdict),4))
    print '\nRegion: {:02d}\n'.format(n_reg), 10*'='
    for lab, ax in zip(labels, axes):
        print '\nEnergy band: {}\n'.format(lab)

        want_err = not(blist_reg and lab in blist_reg)
        if not want_err:
            print 'Blacklisted profile'

        bdict = bdict_fwhm_fit(rdict[lab], fitter, want_err=want_err, **kwargs)
        print 'FWHM = {}, errors {}'.format(bdict['fwhm'], bdict['errs'])

        rdu.plot_fwhm_fit(bdict, ax)
        rdict[lab] = bdict
    plt.show()

    plt.figure(figsize=(8,6))
    for lab in labels:
        bdict = rdict[lab]
        fwhm = bdict['fwhm']
        if np.isfinite(fwhm):
            nerr, perr = bdict['errs']
            en = float(re.split('-', lab)[0])
            plt.errorbar([en], [fwhm], yerr=[[nerr], [perr]], fmt='bs')
    fplot('Energy (keV)', 'FWHM (arcsec)', ax=plt.gca(),
          axargs=[0., 7., plt.ylim()[0], plt.ylim()[1]])
    plt.show()

    return rdict

def bdict_fwhm_fit(bdict, fitter, cap=False, sub_bkg=True, **kwargs):
    """Compute and store best fit FWHM in band dictionary

    Overwrites any previous fit information

    Input:
        bdict: band dictionary
        fitter: fit routine w/ input f(x, y, y_err); output popt, pcov, f, bkg
                (just augments the output of sp.optimize.curve_fit)
        cap (boolean): apply a cap at data max for FWHM computation?
        sub_bkg (boolean): subtract background in FWHM computation?
        Remaining **kwargs (want_err, dx, xconst) sent to ffwhm.get_fwhm_all
    Output:
        new band dictionary.  also modifies input bdict anyways
    """
    x, y, y_err, c = bdict['data']
    popt, pcov, f, bkg = fitter(x[c:], y[c:], y_err[c:])

    if sub_bkg:
        f_fwhm = lambda x, *pars: f(x,*pars) - bkg
        y = y - bkg
    else:
        f_fwhm = f
    y_cap = np.amax(y[c:]) if cap else None  # Must follow sub_bkg check

    w, lims, errs = ffwhm.get_fwhm_all(x[c:], y[c:], y_err[c:], f_fwhm, popt,
                                       y_cap=y_cap, **kwargs)

    if 'want_err' in kwargs and not kwargs['want_err']:
        w = float('NaN')
        lims = (float('NaN'), float('NaN'))

    f_src = inspect.getsource(f)  # Hackily convert function name to f
    bdict['f'] = re.sub(r'.*def\s*\S+\s*\(', 'def f(', f_src, count=1)
    bdict['pars'] = popt
    bdict['fwhm'] = w
    bdict['lims'] = lims
    bdict['errs'] = errs

    return bdict

# Compute spectrum cut locations for given region dict
# Kind of separate/unrelated functionality, but stuck here for now

def get_cuts(rdict):
    """Cut locations for spectra in arcsec
    Ignores any profile bands where FWHM could not be found
    Returns x_min, x_btw, x_max (NOT STORED IN RDICT)
    """
    x_btw = np.nanmin([bdict['lims'][0] for bdict in rdict.values()])
    x_max = np.nanmax([bdict['lims'][1] for bdict in rdict.values()])

    x_min = -1
    for lab, bdict in rdict.items():
        # Only use cuts where FWHM could be fitted
        if np.isfinite(bdict['fwhm']):
            x, y, y_err, c = bdict['data']
            if x_min > x[c] or x_min == -1:
                x_min = x[c]
    if x_min == -1:
        print 'WARNING: could not find fit domain minimum?!'
        x_min = 0

    assert x_min <= x_btw <= x_max

    return x_min, x_btw, x_max

# ==================
# Dump data to files
# ==================

def regdict_dump(rdict, nreg, outroot):
    """Dump regdict to plaintext json and binary"""
    check_dir(outroot)
    with open('{}_{:02d}.pkl'.format(outroot, nreg), 'w') as fobj:
        pickle.dump(rdict, fobj)
    with open('{}_{:02d}.json'.format(outroot, nreg), 'w') as fobj:
        json.dump(rdict, fobj, cls=NumpyJSONEncoder, indent=4)

def speccut_dump(spec_cuts, outroot):
    """Dump spectrum cut information to plaintext and npz"""
    check_dir(outroot)
    fout_dat = '{}_spec_cuts.dat'.format(outroot)
    fout_npz = '{}_spec_cuts.npz'.format(outroot)
    hdr = 'Spectrum cuts, radial pos (arcsec), for processed/fitted profiles'
    np.savetxt(fout_dat, spec_cuts, header=hdr)
    np.savez(fout_npz, cuts=spec_cuts)


if __name__ == '__main__':
    pass
