"""
Perform phabs*powerlaw XSPEC fit to spectra and output
plots and fit parameters
Aaron Tran
2014 June 17
(last modified: 2014 June 17)

Initialize HEASOFT (`heainit`) before running this script!
Run script with 32-bit python (use `arch -i386 python ...`)
Run script from directory of actual spectra (so that XSPEC can resolve relative
links to response and background files.
"""

import argparse
import os
import re

import xspec as xs

def main():
    """Main method"""
    # argparse boilerplate
    parser = argparse.ArgumentParser(description=
             'Apply XSPEC model fits and make plots for set of spectra')
    parser.add_argument('specroot', help='Directory stem for spectra')
    parser.add_argument('plotroot', help='Output stem for plots')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    specroot, pltroot = args.specroot, args.plotroot
    verbose = args.verbose

    # Get number of spectra to update (count grouped spectrum files)
    n = get_nspec(specroot)
    if verbose:
        print '\n{} spectra to process'.format(n)

    # Check plot output directory, create if needed
    check_dir(pltroot)

    # Set up XSPEC
    init_xspec(verbose)

    # Process spectra one by one
    for num in xrange(n):
        # Create paths, with assumptions about filename structure
        # Check that grouped spectrum exists
        grp_path = '{}_src{:d}_grp.pi'.format(specroot, num+1)
        plt_path = '{}_src{:d}_grp.ps'.format(pltroot, num+1)
        log_path = grp_path[:-3] + '.fitlog'
        if not os.path.isfile(grp_path):
            print 'Bad path: {}'.format(grp_path)
            raise Exception('ERROR: file does not exist!')
        if verbose:
            print '\nProcessing file: {}'.format(grp_path)
            print 'Output plot: {}'.format(plt_path)

        process_spectrum_file(grp_path, plt_path, log_path)

    if verbose:
        print '\nDone!'


# ========================
# Methods to control XSPEC
# ========================

def process_spectrum_file(fname, pltname, logname):
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
    
    # Print out a color postscript plot
    xs.Plot.device = pltname + '/cps'
    xs.Plot('ldata', 'residual', 'ratio')

    # Save spectrum/model/fit information to log file
    xs.Xset.logChatter = 10
    s.show()
    model.show()
    xs.Fit.show()
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


# ===========================
# Utility methods for spectra
# ===========================

def check_dir(stem):
    """Check if stem directory exists and create if needed"""
    stemdir = os.path.dirname(stem)
    if not os.path.isdir(stemdir) and stemdir != '':
        print 'stem {} in nonexistent directory'.format(stem)
        print 'Creating directory {}'.format(stemdir)
        os.makedirs(stemdir)


def get_nspec(stem):
    """Number of grouped spectra with given directory stem"""
    pattern = os.path.basename(stem) + r'_src[0-9]+_grp\.pi'
    if os.path.dirname(stem) == '':
        files = [f for f in os.listdir('.') if re.match(pattern,f)]
    else:
        files = [f for f in os.listdir(os.path.dirname(stem))
                 if re.match(pattern, f)]
    return len(files)


if __name__ == '__main__':
    main()

