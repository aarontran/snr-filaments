"""
Code for IPython notebook fitting -- generate and save fits
attempting to abstract a lot of boilerplate stuff out

build_dataf(...), generate_fits(...) work in tandem

factory function f(...) from build_dataf(...) outputs an lmfit.Minimizer()
object (res) w/ lots of information attached to res.params
Relevant properties are:

    res.params.ci           dict of confidence intervals
    res.params.info         dict of best fit parameters w/ keys:
        'conf-intv'             CI value
        'chisqr'                fit chisqr
        'redchi'                fit reduced chisqr
        'ndata'                 res.ndata
        'nvarys'                res.nvarys
        'nfree'                 res.nfree
    res.params.metadata     dict of build_dataf args
        'fit-type'              simp/full
        'fit-kws'               eta2_free, model_kws, etc...
        'err-kws'               manual/lmfit/stderr, eps, adapt, etc...

In the "metadata" pickle:

    fobj_all = {}
    fobj_all['title'] = fobj.title
    fobj_all['kevs'] = fobj.kevs
    fobj_all['data'] = fobj.data
    fobj_all['eps'] = fobj.eps
    fobj_all['inds'] = fobj.inds
    fobj_all['snr-pars'] = fobj.snr.config_log()

    NOTE: kevs, data, eps ALREADY have inds applied
    so do not reapply inds when recreating fitter object.
    Bad design, should be set up to save "original" kevs/data/eps
    then apply inds.  But leave for now...

Aaron Tran
August 2014
"""

from __future__ import division

import cPickle as pickle
from glob import glob
import json
import lmfit  # For type-checking, JSON serialization
import matplotlib.pyplot as plt
import numpy as np
import re

from json_utils import LmfitJSONEncoder
from nice_tables import LatexTable, ListTable
import models
import models_exec as mex
import regdict_utils as rdu
import regparse

from fplot import fplot

# ==============================
# Run fits and save data to disk
# ==============================
def build_dataf(fit_type, conf_intv=0.683, fit_kws=None, err_kws=None):
    """Build function f to perform fits/errs; f(mu, fobj) returns res obj with
    fits, errors, fit properties, etc... (see src)

    f is a wrapper function for fobj.get_fit(...), fobj.get_err(...)
    Please see the methods in `models_exec.py`, to understand how they work
    more deeply.

    Notes: fitting w/ one free param may not give sensible (manual) errors;
    lmfit may balk.  Input options are not all fully tested, but most should
    work (some unit tests would be nice to demonstrate that error finding does
    what I expect it to do).

    Inputs:
        fit_type: 'simp' or 'full'
        conf_intv: please don't set above 0.683 or your life will suck, I guarantee
        fit_kws, *simple and full models*
            eta2_free=True, B0_free=True
            eta2, B0: provide custom (non-grid) initial guesses.
            mu_free=False (I recommend not to change, esp. if using full model!
                           NOT TESTED, likely to hit bugs)
            ab_fit: if provided, will allow damping length ab to vary in fit
                be sure to enable damping in model_kws as well!
                value of ab_fit will override any specified damp_ab value
        fit_kws, *full model only*
            model_kws (dict):
                rminarc, icut, irmax, iradmax, ixmax, irad_adapt, irad_adapt_f,
                idamp, damp_ab, damp_bmin, fgfname, itmax, inmax, irhomax
            scale_covar=False (!) (lmfit.minimize)
            method='leastsq'
            epsfcn, maxfev, factor, diag, ftol, xtol, etc... (scipy leastsq)
        err_kws:
            method='manual' ('lmfit', 'stderr')
            if using 'manual':
                anneal=True, eps=None, adapt=True (passed to mex.one_dir_root)
                sp.optimize.brentq kwargs: xtol, rtol
            if using 'lmfit':
                maxiter=200, prob_func=None (lmfit.conf_interval kwargs)

    Output:
        wrapper function f(mu, fobj) that performs a fit using the supplied
        keywords/settings from this factory function.
        f(...) outputs an lmfit.Minimizer() object (res) w/ information
        attached to res.params
        Relevant properties:
            res.params.ci           dict of confidence intervals
            res.params.info         dict of best fit parameters w/ keys:
                'conf-intv'             CI value
                'chisqr'                fit chisqr
                'redchi'                fit reduced chisqr
                'ndata'                 res.ndata
                'nvarys'                res.nvarys
                'nfree'                 res.nfree
            res.params.metadata     dict of build_dataf args
                'fit-type'              simp/full
                'fit-kws'               eta2_free, model_kws, etc...
                'err-kws'               manual/lmfit/stderr, eps, adapt, etc...
    """
    if fit_kws is None:
        fit_kws = {}
    if err_kws is None:
        err_kws = {}

    def assmb_data(mu, fobj):
        """Assemble fit/error data for a given mu
        Attaches a ton of extra information to res.params object, for later
        retrieval
        """
        # Perform actual fit
        res = fobj.get_fit(mu, fit_type, **fit_kws)

        # Save parameter values and standard errors before error computation
        chisqr, redchi = res.chisqr, res.redchi
        ndata, nvarys, nfree = res.ndata, res.nvarys, res.nfree
        eta2, B0 = res.params['eta2'].value, res.params['B0'].value
        eta2_stderr = res.params['eta2'].stderr
        B0_stderr = res.params['B0'].stderr
        if 'ab' in res.params:
            ab = res.params['ab']
            ab_stderr = res.params['ab'].stderr

        # Compute errors (whether stderr, lmfit, or my manual procedure)
        ci = fobj.get_errs(res, fit_type, ci_vals=(conf_intv,), **err_kws)

        # If eta2/B0 are held fixed, they won't have CI entries
        if 'eta2_free' in fit_kws and not fit_kws['eta2_free']:
            ci['eta2'] = [(conf_intv, eta2), (0., eta2), (conf_intv, eta2)]
        if 'B0_free' in fit_kws and not fit_kws['B0_free']:
            ci['B0'] = [(conf_intv, B0), (0., B0), (conf_intv, B0)]

        # (re)Store all data in res object
        res.chisqr, res.redchi = chisqr, redchi
        res.ndata, res.nvarys, res.nfree = ndata, nvarys, nfree
        res.params['eta2'].value = eta2
        res.params['eta2'].stderr = eta2_stderr
        res.params['B0'].value = B0
        res.params['B0'].stderr = B0_stderr
        if 'ab' in res.params:
            res.params['ab'].value = ab
            res.params['ab'].stderr = ab_stderr

        # Hacky workaround to send fit information downstream
        # Using params because it's pickle-able for IPython's parallel stuff
        # semantically illogical, but it works...
        res.params.ci = ci
        res.params.snr = fobj.snr

        info = {}
        info['conf-intv'] = conf_intv
        info['chisqr'] = res.chisqr
        info['redchi'] = res.redchi
        info['ndata'] = res.ndata
        info['nvarys'] = res.nvarys
        info['nfree'] = res.nfree
        res.params.info = info

        metadata = {}
        metadata['fit-type'] = fit_type
        metadata['fit-kws'] = fit_kws
        metadata['err-kws'] = err_kws
        res.params.metadata = metadata

        return res

    return assmb_data

def generate_fits(f_data, fobj, n_ID, mu_vals, outroot, save=True):
    """Perform fits, save to pkl/json if desired;
    return list of Parameters objects and Fitter()
    """
    p_list = []
    for mu in mu_vals:
        print '{} stdout'.format(fobj.title)
        res = f_data(mu, fobj)
        p_list.append(res.params)
    if save:
        save_fits(p_list, fobj, mu_vals, outroot + '-{:02d}'.format(n_ID))
    return p_list, fobj

def save_fits(p_list, fobj, mu_vals, outroot):
    """Save full model fit parameters to pkl/json
    Use the pkl though, seriously.
    The json is just a backup, for when Python/pickle/etc are obsolute

    This DOES clobber any existing output!

    please see build_dataf for info on stored data!

    Input:
        p_list: list of lmfit.Parameters objects (one per mu value)
        fobj: Fitter() object
        mu_vals: fixed mu values used for fitting
        outroot: output file stem (will suffix stuff to this)
    Output:
        creates pkl of zipped (mu_vals, p_list); each Parameters object in the
        tuple carries ci, info, metadata
        also creates pkl with SNR information
        JSON implementation coming -- have to deconstruct Parameters object
        somehow.
    """
    ci_list = []
    info_list = []
    meta_list = []
    for p, mu in zip(p_list, mu_vals):
        ci_list.append(p.ci)
        info_list.append(p.info)
        meta_list.append(p.metadata)

    # Spew metadata relevant to fits for all mu values
    fobj_all = {}
    fobj_all['title'] = fobj.title
    fobj_all['kevs'] = fobj.kevs
    fobj_all['data'] = fobj.data
    fobj_all['eps'] = fobj.eps
    fobj_all['inds'] = fobj.inds
    fobj_all['snr-pars'] = fobj.snr.config_log()

    # Output everything.  Doesn't check for overwriting
    regparse.check_dir(outroot)

    pkg = [p_list, mu_vals]

    with open('{}-data.pkl'.format(outroot), 'w') as fpkl:
        pickle.dump(pkg, fpkl)
    with open('{}-fobj.pkl'.format(outroot), 'w') as fpkl:
        pickle.dump(fobj_all, fpkl)

    with open('{}-data.json'.format(outroot), 'w') as fjson:
        json.dump([p_list, mu_vals, info_list, ci_list, meta_list],
                  fjson, cls=LmfitJSONEncoder, indent=4)
    with open('{}-fobj.json'.format(outroot), 'w') as fjson:
        json.dump(fobj_all, fjson, cls=LmfitJSONEncoder, indent=4)

# =====================================================
# Load best fits, compute FWHMs or return useful things
# =====================================================

# NOTE code to read/parse saved best fits is not well designed (Oct 28 2014)
# there should be a hierarchy of items passed through / returned
# to allow easy retrievals / operations with Fitter objects

# NOT TOUCHING FOR NOW TO PRESERVE EXISTING FUNCTIONALITY
def load_fit_pkls(inroot, want_fobj=False):
    """Generator that iterates over each fitted region/filament
    Input:
        inroot (str) is base file stem, e.g.,
        ../../data-tycho/fwhms/model-fits/simp-man_err
    Yield, for each region fitted:
        tuple of 1. list of Parameters(), 2. list of mu values, 3., region num
                 4. fobj_all dict w/ fitting data, if desired
    """
    npkls = len(glob('{}-[0-9]*-data.pkl'.format(inroot)))
    for n in xrange(npkls):
        fname = '{}-{:02d}-data.pkl'.format(inroot, n+1)  # Enforce ordering
        ffobj = '{}-{:02d}-fobj.pkl'.format(inroot, n+1)
        with open(fname, 'r') as fpkl:
            p_list, mu_vals = pickle.load(fpkl)
        if want_fobj:
            with open(ffobj, 'r') as fpkl:
                fobj_dict = pickle.load(fpkl)
            yield p_list, mu_vals, n+1, fobj_dict
        else:
            yield p_list, mu_vals, n+1

def load_fit_pkls_new(inroot):
    """Generator iterating over each fitted region/filament
    Input:
        inroot (str) is base file stem, e.g.,
        ../../data-tycho/fwhms/model-fits/simp-man_err
    Yield, for each region fitted:
        list of Parameters()
        list of mu values
        Fitter() object
        region number
    """
    npkls = len(glob('{}-[0-9]*-data.pkl'.format(inroot)))
    for n in xrange(npkls):  # Enforces ordering
        p_list, mu_vals, fobj = load_fit_pkl(inroot, n+1)
        yield p_list, mu_vals, fobj, n+1

def load_fit_pkl(inroot, n):
    """Get best fit information/usables from specified fit pickle

    Input:
        inroot (str): base file stem for best fit pickles
        n (int): region number for inroot
    Output:
        (p_list, mu_list, fobj)

        p_list is list of lmfit.Parameters objects
        mu_vals contains all fitted mu values
        fobj is a regenerated models_exec.Fitter object
            (assumed same for all mu values considered)
            set to verbose by default
    """

    fname = '{}-{:02d}-data.pkl'.format(inroot, n)
    ffobj = '{}-{:02d}-fobj.pkl'.format(inroot, n)

    with open(fname, 'r') as fpkl:
        p_list, mu_vals = pickle.load(fpkl)
    with open(ffobj, 'r') as fpkl:
        fdict = pickle.load(fpkl)

    p0 = p_list[0]  # Fitter object assumed same for all mu values
    fobj = mex.Fitter(p0.snr, fdict['kevs'], fdict['data'], fdict['eps'],
                      None,  # No table needed
                      inds=None, verbose=True)
    fobj.title = fdict['title']
    return p_list, mu_vals, fobj

def best_fit_fwhms(p_list, mu_vals, fobj, kevs_calc=None, mu_calc=None,
                   get_mask=False, verbose=True, **kwargs):
    """Compute model FWHMs from parameter list, mu values, and Fitter() object
    (obtained from best fit pickle methods)

    Input:
        p_list: list of lmfit.Parameters() objects, one per mu in mu_vals
        mu_vals: list of mu values fitted
        fobj: models_exec.Fitter() object used for fits in p_list/mu_vals
        kevs_calc (list): energies (keV) at which to compute model FWHMs
                          if none specified, data keVs used
        mu_calc (list): mu values at which to get best fit parameters
                        if none specified, [1] used
                        if invalid mu values given, they are ignored
        get_mask (bool): return masks used to screen good FWHMs
        verbose (bool): print model kws being used or not
        **kwargs: passed to fobj.width_full(...); overrides stored model kws
    Output:
        kevs_out: list of kevs_calc for each mu, masked where valid FWHMs found
        fwhms_out: list of FWHMs for each mu, masked where valid FWHMs found
        msk_out: list of boolean masks (size matches kevs_calc) for each mu
            (only if get_mask is True)

        if only one mu value is desired (or present),
        kevs, fwhms, mask are returned directly (automatically unpacked...)

        WARNING: if get_prfs=True is passed in, mask will not be applied.
        function will return kevs, model_outputs as below
        but model outputs will be 3 tuple of FWHMs, intensities, r-coord grids
    """
    if kevs_calc is None:
        kevs_calc = fobj.kevs
    if mu_calc is None:
        mu_calc = [1]

    getting_prfs = 'get_prfs' in kwargs and kwargs['get_prfs']
    # Assumes that get_prfs will not be enabled in saved best fit kws

    kevs_out = []
    model_out = []
    masks_out = []

    fobj.set_kevs(kevs_calc)  # Usage of inds in fobj is very sketchy... beware

    for p, mu in zip(p_list, mu_vals):
        if mu not in mu_calc:
            continue

        # Set up and run fit
        try:
            mkws = p.metadata['fit-kws']['model_kws']
        except KeyError:
            mkws = {}
        eta2 = p['eta2'].value
        B0 = p['B0'].value
        if verbose:
            print 'Original fit kwargs:', mkws
        mkws.update(kwargs)
        if verbose:
            print' Using fit kwargs:', mkws

        # width_cont_return may be FWHMs or (FWHMS, intensities, rgrids)
        width_cont_return = fobj.width_full(mu, eta2, B0, **mkws)

        if getting_prfs:  # NO MASKING IS APPLIED
            kevs_out.append(kevs_calc)
            model_out.append(width_cont_return)
        else:
            fwhms_m = width_cont_return
            mask = fwhms_m < 0.99*p.snr.rsarc  # Mask gives pts w/ valid FWHMs

            kevs_out.append(kevs_calc[mask])
            model_out.append(fwhms_m[mask])
            masks_out.append(mask)

    if len(mu_calc) == 1:
        kevs_out = kevs_out[0]
        model_out = model_out[0]
        if not getting_prfs:
            masks_out = masks_out[0]

    if getting_prfs:
        return kevs_out, model_out

    if get_mask:  # Redundant code, but makes control flow clearer...
        return kevs_out, model_out, masks_out
    else:
        return kevs_out, model_out

# =========================
# Functions to display fits
# =========================
def generate_tabs(p_list, title, mu_vals):
    """Parse list of Parameters() objects, mu values, w/ Fitter()
    to generate quick tables for IPython, LaTeX display.
    Input:
        p_list: output of generate_fits, list of Parameters()
        fobj: Fitter() object
        mu_vals: list of mu values, matched to p_list
    Output:
        table, ltab (IPython, LaTeX tables)
    """

    table = ListTable()
    table.append(['mu', 'eta2', 'B0', 'chisqr'])

    latex_hdr = [r'$\mu$ (-)', r'$\eta_2$ (-)', r'$B_0$ ($\mu$G)', r'$\chi^2$']
    latex_cols = ['{:0.2f}', 2, 2, '{:0.4f}']
    ltab = LatexTable(latex_hdr, latex_cols, title, prec=4)

    for p, mu in zip(p_list, mu_vals):
        #print p.ci
        #print p.info
        #print p.metadata
        # Pull in full errors from p.ci  (get_ci_errors order: pos/neg errs)
        p['B0'].fullerr = mex.get_ci_errors(p.ci, 'B0', p.info['conf-intv'])
        p['eta2'].fullerr = mex.get_ci_errors(p.ci, 'eta2',p.info['conf-intv'])

        # Edit B0 parameter for table formatting
        B0_min, p['B0'].min = p['B0'].min, None  # scrub but save B0 limits
        B0_max, p['B0'].max = p['B0'].max, None
        p['B0'].value *= 1e6
        p['B0'].stderr *= 1e6
        p['B0'].fullerr = 1e6 * np.array(p['B0'].fullerr)  # Hackish workaround

        # Build iPython table
        tr = ['{:0.2f}'.format(mu)]
        tr.append(('{0.value:0.3f} +{1[0]:0.2f}/-{1[1]:0.2f} '
                   '(std: &plusmn; {0.stderr:0.3f})').format(
                  p['eta2'], p['eta2'].fullerr))
        tr.append(('{0.value:0.3g} +{1[0]:0.3g}/-{1[1]:0.3g} '
                   '(std: &plusmn; {0.stderr:0.3f})').format(
                  p['B0'], p['B0'].fullerr))
        tr.append('{:0.4f}'.format(p.info['chisqr']))
        table.append(tr)

        # Build LaTeX table
        ltab.add_row(mu, p['eta2'].value, p['eta2'].fullerr[0],
            p['eta2'].fullerr[1], p['B0'].value, p['B0'].fullerr[0],
            p['B0'].fullerr[1], p.info['chisqr'])

        # Reset B0 value/errors, in case list is used elsewhere
        p['B0'].value *= 1e-6
        p['B0'].stderr *= 1e-6
        p['B0'].fullerr *= 1e-6
        p['B0'].min = B0_min  # Must follow value resetting
        p['B0'].max = B0_max

    return table, ltab

def generate_plots(p_list, fobj, mu_vals, fmt_vals, ax=None):
    """Makes plots.  Yeah.  Logical follow-on to generate_tabs"""
    if ax is None:
        plt.figure()
        ax = plt.gca()
    ax.errorbar(fobj.kevs, fobj.data, yerr=fobj.eps, fmt='ok', label='Data')

    for p, mu, fmt in zip(p_list, mu_vals, fmt_vals):

        fit_type = p.metadata['fit-type']
        fit_kws = p.metadata['fit-kws']

        if fit_type == 'simp':
            kevs_m = np.linspace(fobj.kevs[0]-0.2, fobj.kevs[-1]+0.2, 100)
            fwhms_m = models.width_dump(p, kevs_m, fobj.snr)

        elif fit_type == 'full':
            fill_kevs = np.linspace(fobj.kevs[0]-0.2, fobj.kevs[-1]+0.2, 5)
            kevs_m = np.sort(np.hstack((fobj.kevs, fill_kevs)))
            if fit_kws is not None and 'model_kws' in fit_kws:
                model_kws = fit_kws['model_kws']
            else:
                model_kws = {}
            fwhms_m = models.width_cont(p, kevs_m, fobj.snr, verbose=False,
                                        **model_kws)
        ax.plot(kevs_m, fwhms_m, fmt, label=r'$\mu = {:0.2f}$'.format(mu))

    # Prepare composite plot
    fplot('Energy (keV)', 'FWHM (arcsec)', axargs='tight')
    ax.legend(loc='best')
    ax.set_title(fobj.title)
    return ax

# ========================
# Functions to set-up fits
# ========================

def init_data(rroot, kevs, labels):
    """Prepared regdict data for model fitting
    Averages asymmetric FWHM errors
    Input:
        rroot (str): region dict filestem
        kevs (list): energy band values, assoc. w/ labels
        labels (list): energy band labels, strings
    Output:
        dictionary keyed by region number to (kevs, data, eps, inds)
    """
    data_all = {}

    for rdict, n_reg in rdu.regdict_load(rroot):
        data = np.array([rdict[lab]['fwhm'] for lab in labels])
        eps = np.mean([rdict[lab]['errs'] for lab in labels], axis=1)
        inds = (np.where(np.isfinite(data)))[0]
        data_all[n_reg] = kevs, data, eps, inds

    return data_all

# TODO if we ever use averaged FWHMs,
# splice this out into a separate script
# make a BRAND NEW set of region dictionaries + azimuth angles
# then, perform analysis as usual -- instead of treating specially here
def init_data_avg(data_regs, flmts, avg='a'):
    """Average region FWHM data that has already been retrieved

    Currently reporting standard error of mean as errors
    but, for small sample sizes, we must use Student's t-distribution.
    For 70% confidence interval, multiply by 1.963 for dof=1 (n=2)
    For 90% confidence interval, multiply by 6.3 for dof=1 (n=2) (!!!)

    Input:
        data_all (dict): keyed by region number, w/ kevs/data/eps/inds (use init_data)
        flmts (dict): filament number paired w/ constituent region numbers
        avg (str): 'a' or 'g' for arithmetic, geometric averages & errors
    Output:
        new dict, keyed by filament numbers, w/ averaged data
        also includes n_pts (number of points averaged in each band)
    """
    data_avgd = {}

    for n_flmt, n_regs in flmts.items():
        kevs = []
        data_all = []
        for m in n_regs:
            kevs_r, data_r, eps_r, inds_r = data_regs[m]
            if len(kevs) > 0 and not np.array_equal(kevs, kevs_r):
                raise Exception('Discrepancy in stored kevs')
            kevs = kevs_r
            data_all.append(data_r)

        data_all = np.array(data_all)
        n_pts = np.sum(~np.isnan(data_all),axis=0)  # Number of usable FWHMs in each energy band
        inds = np.where(n_pts)[0]

        # Lots of unsafe computation -- but, if n=0 we will throw it out anyways
        if avg == 'a':
            data = np.nanmean(data_all, axis=0)
            std = np.nanstd(data_all, axis=0, ddof=1)
            std[np.where(np.isnan(std))] = np.nanmax(std)  # Catch n=1 cases (stdev ill-defined)
            eps = std / np.sqrt(n_pts)  # Standard error of the mean
        elif avg == 'g':
            # Compute errors in log space first, before setting NaNs in data to 1
            std_log = np.nanstd(np.log(data_all), axis=0, ddof=1)
            std_log[np.where(np.isnan(std_log))] = np.nanmax(std_log)  # Catch n=1 cases (stdev ill-defined)
            eps_log = std_log / np.sqrt(n_pts)  # Convert to standard err in log space

            data_all[np.isnan(data_all)] = 1  # Identity for multiplying
            data = np.power(np.prod(data_all, axis=0), 1./n_pts)

            # Convert back to original parameter space
            eps_upper = np.exp(np.log(data) + eps_log) - data
            eps_lower = data - np.exp(np.log(data) - eps_log)
            # We need to provide a symmetric error for fitting
            eps = np.maximum(eps_upper, eps_lower)

        data_avgd[n_flmt] = kevs, data, eps, inds, n_pts

    return data_avgd


if __name__=='__main__':
    pass

