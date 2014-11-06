"""
Perform XSPEC fits to characterize spectral variations in filaments
Aaron Tran

Base fit is phabs*po with nH = 0.7 frozen.  Then we have some options:
0. phabs*po alone
1. phabs*po with S line excised (cut 2.3-2.6 keV)
2. phabs*po with S line fit

I will also run manual fits to consider other lines...
Those may not get completely logged though
"""

import argparse
import json
import os
import numpy as np

import xspec as xs

import regparse
import spec_fit_utils as xsutils

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
                                         '1=excise S line, '
                                         '2=excise S, Ar lines, '
                                         '3=fit S line'), type=int)
    parser.add_argument('plotroot', help='Output stem for plots')
    parser.add_argument('fitproot', help='Output stem for fit logs, data')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    args = parser.parse_args()
    specroot, pltroot, fitproot = args.specroot, args.plotroot, args.fitproot
    ftype = args.fittype
    verbose = args.verbose

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

        spec = xs.Spectrum(grp_path)
        model = run_fit(spec, ftype)  # XSPEC fitting here

        # Collective dump and clean-up
        xsutils.plot_dump_data(spec, model, plt_path, log_root+'.npz')
        xsutils.dump_fit_log(spec, model, log_root+'.log')
        xsutils.dump_fit_dict(xsutils.fit_dict(spec, model, want_err=False),
                              log_root+'.json')
        xs.AllData.clear()
        xs.AllModels.clear()

    if verbose:
        print '\nDone!'


def run_fit(spec, ftype=0):
    """Fit spectrum to desired XSPEC model (magic numbers galore!)
    Here, twiddle with desired fit guesses, freezing/thawing, etc

    Input
        fname (str): file of single spectrum
        ftype (int): 0,1,2 -- type of fit to perform
            0: fit to absorbed power law
            1: fit to absorbed power law, ignoring 2.3-2.6 keV data
            2: fit to absorbed power law + gaussian at 2.45 keV
    Output
        2-tuple of xs.Spectrum, xs.Model objects
    """

    spec.ignore('**-2.0, 7.0-**')  # Tail only for spectral variation!

    model = xs.Model('phabs*powerlaw')
    model.phabs.nH = 0.7
    model.phabs.nH.frozen = True
    model.powerlaw.PhoIndex = 1
    model.powerlaw.norm = 1
    xs.Fit.perform()  # First round fit to phabs*po alone, regardless of ftype

    if ftype == 1:  # Excise S line

        spec.ignore('2.3-2.6')
        xs.Fit.perform()
        spec.notice('2.3-2.6')

    elif ftype == 2:  # Excise S and Ar lines

        spec.ignore('2.3-2.6, 3.0-3.1')
        xs.Fit.perform()
        spec.notice('2.3-2.6, 3.0-3.1')

    elif ftype == 3:  # Fit S line

        # Get parameter values from current phabs*po model
        plist = [p.values[0] for p in [model(1), model(2), model(3)]]
        plist.extend([2.45, 2e-2, 5e-7]) # LineE, Sigma, norm

        # Make new model w/ Si line, restore old parameters and add guesses
        model = xs.Model('phabs*(powerlaw + gaussian)')
        model.setPars(*plist)
        model.phabs.nH.frozen=True  # Must reissue freeze command

        model.gaussian.Sigma.frozen=True
        model.gaussian.LineE.frozen=True
        xs.Fit.perform()

        model.gaussian.Sigma.frozen=False
        xs.Fit.perform()
        model.gaussian.LineE.frozen=False
        xs.Fit.perform()

    return model


if __name__ == '__main__':
    main()


