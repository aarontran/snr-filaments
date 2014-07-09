"""
Perform phabs*powerlaw XSPEC fit to spectra and output
plots and fit parameters
Aaron Tran
2014 June 17
(last modified: 2014 July 9)

Initialize HEASOFT (`heainit`) before running this script!
Run script with 32-bit python (use `arch -i386 python ...`)
Run script from directory of actual spectra (so that XSPEC can resolve relative
links to response and background files.
"""

import argparse
import json
import os
import re
import numpy as np

import xspec as xs

import regparse

def main():
    """Main method"""
    # argparse setup
    parser = argparse.ArgumentParser(description=
             ('Apply XSPEC model fits and output spectra/fits/data. '
              'Run script from spectra-containing directory, '
              'so XSPEC can locate response/background files.'))
    parser.add_argument('specroot', help='Directory stem for spectra')
    parser.add_argument('fittype', help=('Type of fit to apply; '
                                         '0=phabs*po, '
                                         '1=excise Si line, '
                                         '2=fit Si line'), type=int)
    parser.add_argument('plotroot', help='Output stem for plots')
    parser.add_argument('fitproot', help='Output stem for fit logs, data')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    specroot, pltroot, fitproot = args.specroot, args.plotroot, args.fitproot
    ftype = args.fittype
    verbose = args.verbose

    # Get number of spectra to fit (count grouped spectrum files)
    n = regparse.count_files_regexp(specroot + r'_src[0-9]+_grp\.pi')
    if verbose:
        print '\n{} spectra to process'.format(n)

    # Check plot output directory, create if needed
    regparse.check_dir(pltroot, verbose)
    regparse.check_dir(fitproot, verbose)

    # Set up XSPEC
    init_xspec(verbose)

    for num in xrange(n):
        # Create paths, with assumptions about filename structure
        # Check that grouped spectrum exists
        grp_path = '{}_src{:d}_grp.pi'.format(specroot, num+1)
        plt_path = '{}_src{:d}_grp.ps'.format(pltroot, num+1)
        log_root = '{}_src{:d}_grp'.format(fitproot, num+1)
        if not os.path.isfile(grp_path):
            print 'Bad path: {}'.format(grp_path)
            raise Exception('ERROR: file does not exist!')
        if verbose:
            print '\nProcessing file: {}'.format(grp_path)
            print 'Output plot: {}'.format(plt_path)

        # All XSPEC interaction here
        spec, model = perform_fit(grp_path, ftype)
        output_fit(spec, model, plt_path, log_root, ftype)

    if verbose:
        print '\nDone!'


# ============================
# Subroutines to control XSPEC
# ============================

def output_fit(s, m, pltname, logroot, ftype):
    """Plot and print results of fit.  Also computes errors."""

    # Print color postscript plot
    xs.Plot.device = pltname + '/cps'
    xs.Plot('ldata', 'residual', 'ratio')

    # Save spectrum/model/fit information to log file
    logFile = xs.Xset.openLog(logroot+'.log')
    xs.Xset.logChatter = 10
    s.show()
    m.show()
    xs.Fit.show()
    if ftype == 2:  # Fitting line
        xs.AllModels.eqwidth(3)
    #xs.Fit.error('1-{}'.format(m.nParameters))  # 90% conf. limit errors
    xs.Xset.logChatter = 0

    # Save fit information to JSON file
    fdict = {}
    fdict['fname'] = s.fileName
    fdict['ftype'] = ftype
    fdict['fitstat'] = (xs.Fit.statMethod, xs.Fit.statistic)
    fdict['dof'] = xs.Fit.dof
    if ftype == 2:
        fdict['eqwidth'] = s.eqwidth[0]
    fdict['pars']={}
    for i in xrange(m.nParameters):
        p = m(i+1)  # xs.Parameter object
        pardict = fdict['pars'][p.name] = {}
        pardict['value'] = p.values[0]
        pardict['error'] = p.error
    
    with open(logroot+'.json','w') as fj:
        json.dump(fdict, fj, indent=2)  # Pretty print


    # Save spectrum data and fit to .npz
    np.savez(logroot,
             x = xs.Plot.x(),
             xE = xs.Plot.xErr(),
             y = xs.Plot.y(),
             yE = xs.Plot.yErr(),
             m = xs.Plot.model(),
             bkg = xs.Plot.backgroundVals())

    # Clean up
    xs.AllData.clear()
    xs.AllModels.clear()
    xs.Xset.closeLog()


def perform_fit(fname, ftype=0):
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
    # Import and fit data
    s = xs.Spectrum(fname)
    s.ignore('**-0.5, 7.0-**')  # 0.5 keV to verify no oxygen line

    # Set up model
    model = xs.Model('phabs*powerlaw')
    model.phabs.nH = 0.6
    model.powerlaw.PhoIndex = 1
    model.powerlaw.norm = 1

    # Fit model
    xs.Fit.nIterations = 2000
    xs.Fit.perform()  # First round fit, regardless of ftype

    # Additional steps to deal with silicon line
    if ftype == 1:  # Excise line
        s.ignore('1.7-2.0')
        xs.Fit.perform()
    elif ftype == 2:  # Fit line
        model = xs.Model('phabs*(powerlaw + gaussian)')
        plist = [m.values[0] for m in [model(1), model(2), model(3)]]
        plist.extend([1.85, 5e-2, 5e-6]) # LineE, Sigma, norm
        model.setPars(*plist)
        # Set soft/hard limits on Si line gaussian Sigma/LineE
        model.gaussian.Sigma.values = [5e-2, 1e-4, 1e-4, 2e-3, 0.07, 0.1]
        model.gaussian.LineE.frozen=True
        xs.Fit.perform()  # Fit with frozen LineE
        model.gaussian.LineE.values = [1.85, 0.001, 1.75, 1.8, 1.9, 1.95]
        model.gaussian.LineE.frozen=False
        xs.Fit.perform()  # Now fit with unfrozen LineE

        # Finally, run fit without hard/soft limits
        par4 = model.gaussian.LineE.values[0]
        par5 = model.gaussian.Sigma.values[0]
        model.gaussian.LineE.values = [par4, 0.01*par4, 0., 0., 1e6, 1e6]
        model.gaussian.Sigma.values = [par5, 0.01*par5, 0., 0., 10., 20.]
        xs.Fit.perform()

    return s, model


def init_xspec(verbose=False):
    """Set global parameters for XSPEC plots, fits
    Returns XSPEC model for fitting
    """
    xs.Plot.xAxis = 'keV'
    xs.Plot.yLog = True
    xs.Plot.background = True

    xs.Xset.addModelString('neivers', '2.0')
    if verbose:
        xs.Xset.chatter = 10
    else:
        xs.Xset.chatter = 10  # hah.
    xs.Xset.logChatter = 0


if __name__ == '__main__':
    main()

