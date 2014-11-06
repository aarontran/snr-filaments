"""
Perform various XSPEC fits and output plots and fit parameters
Aaron Tran

Initialize HEASOFT (`heainit`) before running this script!
Run script from directory of actual spectra (so that XSPEC can resolve relative
links to response and background files).

NOTE: current starting guesses are generally tuned for Tycho, may need to be
adjusted for other SNRs
"""

import argparse
import json
import os
import numpy as np

import xspec as xs

import regparse
import spec_fit_utils as xsutils

SRCUTLOG_PATH_DEFAULT = '/Users/atran3/snr-research/code/srcutlog/'
ALPHA_DEFAULT = 0.58

def main():
    """Parse input and tell XSPEC what to do"""
    parser = argparse.ArgumentParser(description=
             ('Apply XSPEC model fits to all spectra w/ given filestem '
              'and output best fit plots, parameters. '
              'Run from spectra-containing directory '
              'so XSPEC can locate response/background files.'))
    parser.add_argument('specroot', help='Directory stem for spectra')
    parser.add_argument('fittype', help=('Type of fit to apply: '
                                         '0=phabs*po, '
                                         '1=excise Si,S lines, '
                                         '2=fit Si,S lines, '
                                         '3=phabs*srcutlog'), type=int)
    parser.add_argument('plotroot', help='Output stem for plots')
    parser.add_argument('fitproot', help='Output stem for fit logs, data')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    parser.add_argument('-e', '--error', action='store_true',
                        help=('Compute/output errors from 90pct confidence lims'
                              ' (may need user input / take extra time).'))
    parser.add_argument('-p', '--path',
                        help=('Override path to srcutlog local model; default '
                              'is: ' + SRCUTLOG_PATH_DEFAULT))
    parser.add_argument('-a', '--alpha',
                        help='alpha value for SNR (default 0.58 for Tycho)')

    args = parser.parse_args()
    specroot, pltroot, fitproot = args.specroot, args.plotroot, args.fitproot
    ftype = args.fittype
    verbose = args.verbose
    wanterr = args.error
    srcutlog_path, alpha = args.path, args.alpha

    # Check and create new directories if needed
    regparse.check_dir(pltroot, verbose)
    regparse.check_dir(fitproot, verbose)

    n = regparse.count_files_regexp(specroot + r'_src[0-9]+_grp\.pi')
    if verbose:
        print '\n{} spectra to process'.format(n)
    if verbose and wanterr:
        print 'Computing errors from 90% confidence limits'

    # Set up XSPEC for fitting
    xsutils.init_xspec(verbose)
    xs.Fit.nIterations = 2000
    if ftype == 3:  # srcutlog fits need to load local model
        if srcutlog_path is None:
            srcutlog_path = SRCUTLOG_PATH_DEFAULT
        if alpha is None:
            alpha = ALPHA_DEFAULT
        if verbose:
            print 'Using alpha = {}'.format(alpha)
        xs.AllModels.lmod('neipkg', srcutlog_path)

    # Fit spectra and dump output on each iteration
    for num in xrange(n):
        grp_path = '{}_src{:d}_grp.pi'.format(specroot, num+1)
        plt_path = '{}_src{:d}_grp.ps'.format(pltroot, num+1)
        log_root = '{}_src{:d}_grp'.format(fitproot, num+1)
        if not os.path.isfile(grp_path):
            print 'Bad path: {}'.format(grp_path)
            raise Exception('ERROR: spectrum does not exist!')
        if verbose:
            print '\nProcessing file: {}'.format(grp_path)
            print 'Output plot: {}'.format(plt_path)

        # Most XSPEC interaction here - the rest is set-up
        s = xs.Spectrum(grp_path)

        # Default "full-spectrum" fits
        s.ignore('**-0.5, 7.0-**')  # 0.5 keV to verify no oxygen line

        if ftype < 3:
            model = phabspo_fit(s, ftype)
        elif ftype == 3:
            model = srcutlog_fit(s, alpha)

        dump_fit(s, model, plt_path, log_root, ftype, wanterr)

    if verbose:
        print '\nDone!'


# ============================
# Subroutines to control XSPEC
# ============================

def dump_fit(spec, model, pltname, logroot, ftype, wanterr):
    """Plot, print, save fitting results.  Computes errors if requested."""

    # Generates plot and saves spectrum, background, model data/errs
    xsutils.plot_dump_data(spec, model, pltname, logroot+'.npz')

    # Must run logging before saving fit information
    # so we can compute eqwidth, errors (and have that logged)
    def run_extras():
        """XSPEC commands that should be logged"""
        if wanterr:  # Get 90% conf. intv errors for all parameters
            xs.Fit.error('1-{}'.format(model.nParameters))
        if ftype == 2:  # fitting line and need eqwidth
            xs.AllModels.eqwidth(3)
            eqwidth_si = spec.eqwidth[0]
            xs.AllModels.eqwidth(4)
            eqwidth_s = spec.eqwidth[0]
            return eqwidth_si, eqwidth_s  # eqwidths not saved elsewhere

    # Generate log file
    extras = xsutils.dump_fit_log(spec, model, logroot+'.log', run_extras)

    # Generate and customize fit information to save
    fdict = xsutils.fit_dict(spec, model, want_err=wanterr)
    if ftype == 2:
        fdict['eqwidth-si'], fdict['eqwidth-s'] = extras
    xsutils.dump_fit_dict(fdict, logroot+'.json')

    # Clean up for next fit
    xs.AllData.clear()
    xs.AllModels.clear()


def srcutlog_fit(s, alpha):
    """Fit spectrum s to srcutlog model.  First fits break/norm with nH frozen
    at 0.6, then lets nH run free as well.  alpha frozen to user specified value
    srcutlog must be already loaded!
    """

    model = xs.Model('phabs*srcutlog')
    model.phabs.nH = 0.6
    model.srcutlog.alpha = alpha
    model.srcutlog.alpha.frozen = True  # Keep frozen in all fits

    par_break = model.srcutlog.__getattribute__('break')  # Yeah...

    model.phabs.nH.frozen = True
    par_break.frozen = False
    model.srcutlog.norm.frozen = False
    xs.Fit.perform()  # Only break, norm free

    model.phabs.nH.frozen = False
    xs.Fit.perform()  # break, norm, nH all free

    return model


def phabspo_fit(s, ftype=0):
    """Fit spectrum to desired XSPEC model (magic numbers galore!)
    Here, twiddle with desired fit guesses, freezing/thawing, etc

    Input
        fname (str): file of single spectrum
        ftype (int): 0,1,2 -- type of fit to perform
            0: fit to absorbed power law
            1: fit to absorbed power law, ignoring 1.7-2 keV data
            2: fit to absorbed power law + gaussian near 1.85 keV
    Output
        2-tuple of xs.Spectrum, xs.Model objects
    """
    # Set up model
    model = xs.Model('phabs*powerlaw')
    model.phabs.nH = 0.6
    model.powerlaw.PhoIndex = 1
    model.powerlaw.norm = 1

    # Fit model
    xs.Fit.perform()  # First round fit to phabs*po alone, regardless of ftype

    # Additional steps to deal with silicon and sulfur lines
    if ftype == 1:  # Excise lines

        s.ignore('1.7-2.0')
        s.ignore('2.3-2.6')
        xs.Fit.perform()
        s.notice('1.7-2.0')  # So they show up in plots
        s.notice('2.3-2.6')

    elif ftype == 2:  # Fit lines

        # Get parameter values from current phabs*po model
        plist = [p.values[0] for p in [model(1), model(2), model(3)]]
        plist.extend([1.85, 2e-2, 5e-6]) # LineE, Sigma, norm

        # Make new model w/ Si line, restore old parameters and add guesses
        model = xs.Model('phabs*(powerlaw + gaussian)')
        model.setPars(*plist)

        model.gaussian.Sigma.frozen=True
        model.gaussian.LineE.frozen=True
        xs.Fit.perform()

        model.gaussian.Sigma.frozen=False
        xs.Fit.perform()
        model.gaussian.LineE.frozen=False
        xs.Fit.perform()

        # Add sulfur line (same hacky parameter porting)
        plist = [p.values[0] for p in [model(1), model(2), model(3), model(4),
                                       model(5), model(6)]]  # Parameters
        plist.extend([2.45, 2e-2, 5e-7]) # LineE, Sigma, norm
        model = xs.Model('phabs*(powerlaw + gaussian + gaussian)')
        model.setPars(*plist)

        model.gaussian_4.Sigma.frozen=True
        model.gaussian_4.LineE.frozen=True
        xs.Fit.perform()

        model.gaussian_4.Sigma.frozen=False
        xs.Fit.perform()
        model.gaussian_4.LineE.frozen=False
        xs.Fit.perform()

    return model


if __name__ == '__main__':
    main()

