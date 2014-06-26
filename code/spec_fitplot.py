"""
Perform phabs*powerlaw XSPEC fit to spectra and output
plots and fit parameters
Aaron Tran
2014 June 17
(last modified: 2014 June 23)

Initialize HEASOFT (`heainit`) before running this script!
Run script with 32-bit python (use `arch -i386 python ...`)
Run script from directory of actual spectra (so that XSPEC can resolve relative
links to response and background files.
"""

import argparse
import os
import re

import xspec as xs

import regparse

def main():
    """Main method"""
    # argparse boilerplate
    parser = argparse.ArgumentParser(description=
             ('Apply XSPEC model fits and make plots for set of spectra. '
              'Script must be run from spectra-containing directory, '
              'so that XSPEC can find response/background files.'))
    parser.add_argument('specroot', help='Directory stem for spectra')
    parser.add_argument('fittype', help=('Type of fit to apply; '
                                         '0=phabs*po, '
                                         '1=excise Si line, '
                                         '2=fit Si line'), type=int)
    parser.add_argument('plotroot', help='Output stem for plots')
    parser.add_argument('fitproot', help='Output stem for fit logs')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    specroot, pltroot, fitproot = args.specroot, args.plotroot, args.fitproot
    ftype = args.fittype
    verbose = args.verbose

    # Get number of spectra to update (count grouped spectrum files)
    n = regparse.count_files_regexp(specroot + r'_src[0-9]+_grp\.pi')
    if verbose:
        print '\n{} spectra to process'.format(n)

    # Check plot output directory, create if needed
    regparse.check_dir(pltroot, verbose)
    regparse.check_dir(fitproot, verbose)

    # Set up XSPEC
    init_xspec(verbose)

    # Process spectra one by one
    for num in xrange(n):
        # Create paths, with assumptions about filename structure
        # Check that grouped spectrum exists
        grp_path = '{}_src{:d}_grp.pi'.format(specroot, num+1)
        plt_path = '{}_src{:d}_grp.ps'.format(pltroot, num+1)
        log_path = '{}_src{:d}_grp.fitlog'.format(fitproot, num+1)
        if not os.path.isfile(grp_path):
            print 'Bad path: {}'.format(grp_path)
            raise Exception('ERROR: file does not exist!')
        if verbose:
            print '\nProcessing file: {}'.format(grp_path)
            print 'Output plot: {}'.format(plt_path)

        process_spectrum_file(grp_path, plt_path, log_path, ftype)

    if verbose:
        print '\nDone!'


# ========================
# Methods to control XSPEC
# ========================

def process_spectrum_file(fname, pltname, logname, ftype=0):
    """Applies desired processing steps to specified spectrum and XSPEC model

    Input
        fname (str): file of single spectrum
    Output
        None
    Side effects
        Saves plot of spectrum, model, and residuals
        Saves fit parameters to plaintext
    """
    # Start logging
    logFile = xs.Xset.openLog(logname)

    # Import and fit data
    s = xs.Spectrum(fname)
    s.ignore('**-0.6, 7.0-**')

    # Set up model
    model = xs.Model('phabs*powerlaw')
    model.phabs.nH = 0.6
    model.powerlaw.PhoIndex = 1
    model.powerlaw.norm = 1

    # Fit model
    xs.Fit.nIterations = 200
    xs.Fit.perform()

    # Deal with the silicon line if desired
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

    
    # Print out a color postscript plot
    xs.Plot.device = pltname + '/cps'
    xs.Plot('ldata', 'residual', 'ratio')

    # Save spectrum/model/fit information to log file
    xs.Xset.logChatter = 10
    s.show()
    model.show()
    xs.Fit.show()
    if ftype == 2:  # Fitting line
        xs.AllModels.eqwidth(3)
    xs.Xset.logChatter = 0

    # Clean up
    xs.AllData.clear()
    xs.AllModels.clear()
    xs.Xset.closeLog()


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
        xs.Xset.chatter = 5
    xs.Xset.logChatter = 0


if __name__ == '__main__':
    main()

