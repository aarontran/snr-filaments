"""
General utilities/tools for running PyXSPEC fits on SNR remnants (this project)
Aaron Tran
November 2014

Useful methods to plot/dump information about spectrum, model, fit, etc
Also to initialize XSPEC etc...

Initialize HEASOFT (`heainit`) before running/importing
Scripts using PyXSPEC need to be run from same directory as spectra (so that
XSPEC can resolve relative links to response/background files)
"""

import json  # json, np for saving fits
import numpy as np

import xspec as xs


def plot_dump_data(spec, model, plt_fname, npz_fname):
    """spectrum, model should be already set-up/fitted, just like in XSPEC"""

    xs.Plot.device = plt_fname + '/cps'  # Save plot to pltname
    xs.Plot('ldata', 'residual', 'ratio')  # A good default

    np.savez(npz_fname,
             x = xs.Plot.x(), xE = xs.Plot.xErr(),
             y = xs.Plot.y(), yE = xs.Plot.yErr(),
             m = xs.Plot.model(),
             bkg = xs.Plot.backgroundVals())

    return None

def dump_fit_dict(fit_dict, json_fname):
    """Just a default method to dump fit dictionary
    Feed output of fit_dict(...) into here, after you have
    customized/manipulated it as desired
    """
    with open(json_fname, 'w') as fobj:
        json.dump(fit_dict, fobj, indent=2)  # Pretty print


def fit_dict(spec, model, want_err=False):
    """Return dictionary of useful fit parameters
    Just factoring out a common task"""

    fdict = {}
    fdict['fname'] = spec.fileName
    fdict['fitstat'] = (xs.Fit.statMethod, xs.Fit.statistic)
    fdict['dof'] = xs.Fit.dof

    # Crappy way to extract component values and errors
    # Hierarchy of keys: 'comps', 'gaussian', 'sigma', 'value'
    fdict['comps']={}
    for cname in model.componentNames:
        comp = eval('model.'+cname)  # Component object, e.g. m.gaussian1

        compdict = fdict['comps'][cname] = {}

        for pname in comp.parameterNames:
            if pname == 'break':  # Hacky workaround for srcutlog break
                p = comp.__getattribute__('break')
            else:
                p = eval('model.'+cname+'.'+pname)  # Parameter object

            compdict[pname] = {}
            compdict[pname]['value'] = p.values[0]
            if want_err:
                compdict[pname]['error'] = p.error

    return fdict


def dump_fit_log(spec, model, log_fname, func_log=None):
    """Initialize XSPEC logging, dumping current spectrum/model/fit info
    If func_log is given, also run that function (assumes no arguments)
    so that the function will trigger log output
    Then close the log
    Input:
        spec, model are XSPEC spectrum, model objects
        log_fname is filename to save plaintext logfile
        func_log is a function taking no arguments, that triggers XSPEC output
    Output:
        output from func_log if specified, else None
        (file written to log_fname)
    """
    logFile = xs.Xset.openLog(log_fname)
    xs.Xset.logChatter = 10

    spec.show()
    model.show()
    xs.Fit.show()

    if func_log is not None:
        ret_val = func_log()
    else:
        ret_val = None

    xs.Xset.logChatter = 0
    xs.Xset.closeLog()

    return ret_val


def init_xspec(verbose=False):
    """Set (default-ish) global parameters for XSPEC plots, fits"""
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
    pass
