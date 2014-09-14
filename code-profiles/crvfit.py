"""
XSPEC-like functionality to freeze fit parameters
Basically extends scipy.optimize.curve_fit

Next time, use lmfit instead -- much easier

Aaron Tran
July 2014 (moved to .py Sept 2014)
"""

import numpy as np
import scipy as sp

def freeze_func(f, pars, ls_freeze):
    """Make frozen function
    Ex: if ls_freeze = [1,2,4], this builds expression:
    lambda x,p0,p3,p5,p6,p7 : f(x,p0,pars[1],pars[2],p3,pars[4],p5,p6,p7)
    """
    lambd_args = 'x,'
    lambd_body = 'f(x,'
    for i in xrange(len(pars)):
        if i in ls_freeze:
            lambd_body += 'pars[{:d}],'.format(i)
        else:
            lambd_args += 'p{:d},'.format(i)
            lambd_body += 'p{:d},'.format(i)

    lambd_args = lambd_args[:-1]
    lambd_body = lambd_body[:-1] + ')'

    # To explain the `f=f, pars=pars`, see http://stackoverflow.com/q/7349785
    # More than one way to skin a cat, though
    lambd_str = 'lambda ' + lambd_args + ', f=f,pars=pars : ' + lambd_body
    return eval(lambd_str)

def run_frz_fit(x, y, y_err, f, pars, ls_freeze, **kwargs):
    """Perform fit with frozen parameters
    func has form func(x, *pars)
    pars is a list of initial guesses
    ls_freeze is a list of parameter indices that are to be frozen.
    """

    if len(ls_freeze) > 0:
        fit_func = freeze_func(f, pars, ls_freeze)
        p_fit = [pars[i] for i in xrange(len(pars)) if i not in ls_freeze]
    else:
        fit_func = f
        p_fit = pars

    popt, pcov = run_fit(x, y, y_err, fit_func, p_fit, **kwargs)

    popt_aug = popt
    if len(ls_freeze) > 0:
        for i in ls_freeze:
            popt_aug = np.insert(popt_aug, i, pars[i])

    return popt, pcov, popt_aug

def run_fit(x, y, y_err, f, pars, **kwargs):
    """Convenience wrapper for sp.optimize.curve_fit"""
    return sp.optimize.curve_fit(f, x, y, p0=pars, sigma=y_err, absolute_sigma=True, **kwargs)

def print_fit_info(x, y, y_err, f, popt, pcov):
    """Pretty-prints fitting output"""
    np.set_printoptions(precision=2)
    print 'Chi-squared: {} with {} degrees of freedom'.format(chi2(x, y, y_err, f, popt), len(y)-len(popt))
    print 'Reduced chi-squared: {}'.format(chi2red(x, y, y_err, f, popt))
    print 'Errors: {}'.format(np.sqrt(np.diag(pcov)))
    #print 'Covariance matrix:', '\n', pcov
    np.set_printoptions(precision=8)  # Restore to default

def chi2(x, y, y_err, f, pars):
    """Chi-squared for a scalar function fitted to data with y-errors

    Inputs:
        x, y, y_err (iterable): the usual suspects
        f (function): function with call signature f(x, *pars)
        pars (list):  list of parameters for function f
    Output:
        chi-squared value (float)
    """
    return np.sum( ((y - f(x, *pars))/y_err)**2 )

def chi2red(x, y, y_err, f, pars):
    """Reduced chi-squared for a scalar function fitted to data with y-errors

    N.B. for a non-linear model fit, assigning 1 degree of freedom
    to each parameter may be an OVERestimate?
    (there may be other caveats/issues at hand)

    Inputs:
        x, y, y_err (iterable): the usual suspects
        f (function): function with call signature f(x, *pars)
        pars (list):  list of parameters for function f
    Output:
        reduced chi-squared value (float)
    """
    return chi2(x, y, y_err, f, pars) / (len(y) - len(pars))


if __name__ == '__main__':
    pass
