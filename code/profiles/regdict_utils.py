"""
Functions to get / count region dictionaries w/ a given file stem.
Not intended to modify or populate them; please see regdict.py for that

This is for code in other folders that needs access to the region dictionaries'
stored FWHMs.  Python imports don't play nice with symlinks.

To avoid circular dependency, don't import regdict.py ...
But, this does need fsmooth to be able to plot the smoothing function.
sigh.  Deal with this later.

Aaron Tran
September 2014
"""

import cPickle as pickle
from glob import glob
#from matplotlib import
import numpy as np

from fplot import fplot
import fsmooth

# =====================
# Working with regdicts
# =====================

# Following plot settings are twiddled for "publication"-quality figures
# May need to rasterize some data, depending on output file sizes

def plot_fwhm_fit(bdict, ax, allfits=True):
    """Plot cut/uncut data with FWHM fit and bounds for given band dict
    if flag allfits=False, don't plot attempted best fit if FWHM is not
    obtained or blacklisted
    """
    x, y, y_err, c = bdict['data']
    wlen, wtype = bdict['wlen'], bdict['wtype']
    f = get_ffit(bdict)
    pars = bdict['pars']
    w = bdict['fwhm']
    lims = bdict['lims']


    # Plot best fit profile
    if allfits or np.isfinite(w):
        x_m = np.linspace(x[c], x[-1], num=100)
        ax.plot(x_m, f(x_m, *pars), '-b', alpha=1)
    # Plot FWHMs if FWHM is resolved
    if np.isfinite(w):
        ax.axvline(lims[0], ls='-', c='k', lw=1, alpha=0.5)
        ax.axvline(lims[1], ls='-', c='k', lw=1, alpha=0.5)
    # Plot data ABOVE best fit
    plot_smoothing_cut(bdict, ax, show_smth=False)

    return ax

def plot_smoothing_cut(bdict, ax, show_smth=True):
    """Plot cut/uncut data, smoothing used to determine cut (optional)
    Adds axis labels and sets limits
    """
    x, y, y_err, c = bdict['data']
    wlen, wtype = bdict['wlen'], bdict['wtype']

    # Plot data outside, inside fit domain
    ax.errorbar(x[:c], y[:c], yerr=y_err[:c], fmt='r.', capsize=0, alpha=0.5,
                markersize=3)
    ax.errorbar(x[c:], y[c:], yerr=y_err[c:], fmt='k.', capsize=0, alpha=1,
                markersize=3)
    if show_smth:  # Smoothing used to determine fit domain
        ax.plot(x, fsmooth.std_smooth(y, wlen, wtype), '-k', alpha = 0.3)

    # IF necessary, use axargs='tight' and ax.set_ylim(0., ax.get_ylim()[1])
    fplot('Radial position (arcsec.)', 'Intensity (a.u.)', ax=ax,
          axargs='tight')
    ax.set_ylim(0., ax.get_ylim()[1])

    return ax

def get_ffit(bdict):
    """Get fit function used for energy band without polluting namespace"""
    exec(bdict['f'])  # This is an implementation detail
    return f

# =========================================
# Functions to load regdicts, spectral cuts
# =========================================

def regdict_load(rdict_root):
    """Glob and generate region dicts to manipulate"""
    for n in xrange(get_nregs(rdict_root)):
        fname = '{}_{:02d}.pkl'.format(rdict_root, n+1)  # Enforce ordering
        with open(fname, 'r') as fobj:
            rdict = pickle.load(fobj)
        yield rdict, n+1

def speccut_load(rdict_root):
    """Read in npz file of speccut data -- array of 3 x positions, arcsec"""
    fname = rdict_root + '_spec_cuts.npz'
    return np.load(fname)['cuts']

def get_nregs(rdict_root):
    """Get number of regions associated w/ rdict_root"""
    return len(glob('{}_*.pkl'.format(rdict_root)))

if __name__ == '__main__':
    pass
