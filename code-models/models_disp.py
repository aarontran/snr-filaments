"""
Code for rich iPython notebook display
and to abstract away a lot of boilerplate kinda stuff

Storage for mega methods + table, LaTeX formatting material
(may be moved/refactored later)

build_dataf(...), generate_tabs_plot(...) work in tandem
classes ListTable, LatexTable for rich output
(ListTable, in particular, is meant to be displayed in an iPython notebook)

Aaron Tran
August 2014
"""


from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import re

import models
import models_exec as mex

from fplot import fplot

# ==================================
# Functions to control, display fits
# ==================================

def build_dataf(fit_type, conf_intv=0.683, fit_kws=None, err_kws=None):
    """Build function f to perform fits/errs; f(mu) returns res obj with
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
            eta2, B0: provide custom (non-grid) initial guesses. For full
                      model, MUST provide guesses for both, else ignored!
            mu_free=False (I recommend not to change, esp. if using full model!
                           likely to hit bugs, not tested)
        fit_kws, *full model only*
            model_kws (dict):
                rminarc, icut, irmax, iradmax, ixmax, irad_adapt, irad_adapt_f
            scale_covar=False (!) (lmfit.minimize)
            method='leastsq'
            epsfcn, maxfev, factor, diag, ftol, xtol, etc... (scipy leastsq)
        err_kws:
            method='manual' ('lmfit', 'stderr')
            if using 'manual':
                anneal=True, eps=None, adapt=True (passed to mex.one_dir_root)
            if using 'lmfit':
                maxiter=200, prob_func=None (lmfit.conf_interval kwargs)
    """
    if fit_kws is None:
        fit_kws = {}
    if err_kws is None:
        err_kws = {}

    def assmb_data(mu, fobj):
        """Assemble fit/error data for a given mu"""
        # Perform actual fit
        res = fobj.get_fit(mu, fit_type, **fit_kws)

        # Save parameter values and standard errors before error computation
        chisqr = res.chisqr
        eta2, B0 = res.params['eta2'].value, res.params['B0'].value
        eta2_stderr = res.params['eta2'].stderr
        B0_stderr = res.params['B0'].stderr

        # Get better errors
        ci_man = fobj.get_errs(res, fit_type, ci_vals=(conf_intv,), **err_kws)

        # Quick patch -- since assmb_data / etc rely on manipulating B0, eta2
        # specifically already
        if 'eta2_free' in fit_kws and not fit_kws['eta2_free']:
            eta2_err = (0., 0.)
        else:
            eta2_err = mex.get_ci_errors(ci_man, 'eta2', ci_val=conf_intv)

        if 'B0_free' in fit_kws and not fit_kws['B0_free']:
            B0_err = (0., 0.)
        else:
            B0_err = mex.get_ci_errors(ci_man, 'B0', ci_val=conf_intv)

        # (re)Store all data in res object
        res.chisqr = chisqr
        res.params['eta2'].value = eta2
        res.params['eta2'].stderr = eta2_stderr
        res.params['eta2'].fullerr = eta2_err
        res.params['B0'].value = B0
        res.params['B0'].stderr = B0_stderr
        res.params['B0'].fullerr = B0_err

        # Hacky workaround to send fit information downstream
        # Using params because it's pickle-able for iPython's parallel stuff
        # semantically illogical, but it works...
        res.params.fit_type = fit_type
        res.params.fit_kws = fit_kws

        return res

    return assmb_data

def generate_tabs(f_data, fobj, title, mu_vals):
    """f_data takes mu, Fitter(...) as arguments, return lmfit.Minimizer() w/fits+errors
    Yes it's stupid but it will do the job...
    LaTeX tables will NOT include standard errors
    """
    table = ListTable()
    table.append(['mu', 'eta2', 'B0', 'chisqr'])
    ltab = LatexTable([r'$\mu$ (-)', r'$\eta_2$ (-)', r'$B_0$ ($\mu$G)',
                       r'$\chi^2$'],
                      ['{:0.2f}', 2, 2, '{:0.4f}'], title)
    plist = []

    for mu in mu_vals:
        # Get fit/error data
        print '{} computation stdout'.format(title)
        res = f_data(mu, fobj)
        p = res.params

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
        tr.append('{:0.4f}'.format(res.chisqr))
        table.append(tr)

        # Build LaTeX table
        ltab.add_row(mu,
                     p['eta2'].value, p['eta2'].fullerr[0], p['eta2'].fullerr[1],
                     p['B0'].value, p['B0'].fullerr[0], p['B0'].fullerr[1],
                     res.chisqr)

        # Reset B0 value, errors
        p['B0'].value *= 1e-6
        p['B0'].stderr *= 1e-6
        p['B0'].fullerr *= 1e-6
        p['B0'].min = B0_min  # Must follow value resetting
        p['B0'].max = B0_max
        plist.append(p)

    return table, ltab, plist, fobj  # plist, fobj for convenience...

def generate_plots(plist, fobj, title, mu_vals, fmt_vals):
    """Makes plots.  Yeah.  Logical follow-on to generate_tabs"""
    plt.figure()
    ax = plt.gca()
    ax.errorbar(fobj.kevs, fobj.data, yerr=fobj.eps, fmt='ok',
                label='Data')

    for p, mu, fmt in zip(plist, mu_vals, fmt_vals):
        if p.fit_type == 'simp':
            kevs_m = np.linspace(fobj.kevs[0]-0.2, fobj.kevs[-1]+0.2, 100)
            fwhms_m = models.width_dump(p, kevs_m, fobj.snr)
        elif p.fit_type == 'full':
            fill_kevs = np.linspace(fobj.kevs[0]-0.2, fobj.kevs[-1]+0.2, 5)
            kevs_m = np.sort(np.hstack((fobj.kevs, fill_kevs)))
            if p.fit_kws is not None and 'model_kws' in p.fit_kws:
                model_kws = p.fit_kws['model_kws']
            else:
                model_kws = {}
            fwhms_m = models.width_cont(p, kevs_m, fobj.snr, verbose=False,
                                        **model_kws)
        ax.plot(kevs_m, fwhms_m, fmt, label=r'$\mu = {:0.2f}$'.format(mu))

    # Prepare composite plot
    fplot('Energy (keV)', 'FWHM (arcsec)', axargs='tight')
    ax.legend(loc='best')
    ax.set_title(title)
    plt.show()

# ===============
# Utility classes
# ===============

class ListTable(list):
    """ Overridden list class which takes a 2-dimensional list of
    the form [[1,2,3],[4,5,6]], and renders an HTML Table in ipynb.
    Source: http://calebmadrigal.com/display-list-as-table-in-ipython-notebook/
    """

    def _repr_html_(self):
        html = ["<table style='font-family: monospace; font-size: 10pt'>"]
        for row in self:
            html.append("<tr>")

            for col in row:
                html.append("<td>{0}</td>".format(col))

            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)


class LatexTable(object):
    """Format numpy array data with errors in a nice LaTeX'd table.
    Designed for LaTeX package 'booktabs' (what about aastex deluxetable?)

    This duplicates a LOT of functionality from `matrix2latex.py`.
    Possible implementation (using doctest syntax):
    >>> from matrix2latex import matrix2latex as m2ltx
    >>> algnmt = '@{}rr@{ $\pm$ }lr@{ $\pm$ }lr@{}'
    >>> hr = '''[[''] + 5 * [title],
    ...       [r'$\mu$ (-)'] + 2 * [r'$\eta_2$ (-)] + 2 * [r'$B_0$ ($\mu$G)']
    ...        + [r'$\chi^2$']]'''
    >>> fmtcol = 'stuff'  # Here parse sigfigs, depending on error
    >>> # note: fmtcol varies w/ row, not implemented in m2ltx
    >>> m2ltx(data, headerRow = hr, formatColumn = fmtcol, alignment = algnmt)

    I'd like to extend matrix2latex by adding support for data w/errors;
    * use errors to compute # of sigfigs and supply that as an argument
    * use \e in LaTeX to change exponential notation formatting, rather than
      manually editing floats

    For this class, some to-dos (low priority)
    1. customize number of columns, +/- formatting, number formatting
    2. allow concatenating tables together (horizontally)
    3. customize header/footer/layout control

    Options for columns:
    0: plain number, default {:0.3g} formatting
    1: number w/ +/- error, uses +/- symbol as column separator to align
    2: number w/ distinct +/- errors, aligns to left
    str: any string format specifier you like (single right aligned column)

    add_row then accepts an arbitrary number of columns, which should
    correspond to the length of your input column "types" array, +1 for each
    "type 1" column you have and +2 for each "type 2" column
    """

    # Initializes: header, labels, rows, footer; rspec, types
    def __init__(self, labs, types, title):
        """Types specify column types; 0,1,2 = no/single/double +/- errors"""

        self.types = types

        # Build table/column spec
        hlst = [r'\begin{tabular}{@{}']
        for t in types:
            if t == 1:
                hlst.append(r'r@{ $\pm$ }l')  # Use column sep for alignment
            elif t == 2:
                hlst.append('l')  # Align output to values, not errors
            else:
                hlst.append('r')
        hlst.append(r'@{}}')
        self.header = [''.join(hlst)]

        # Table header -- top header, labels, midrule
        self.header.append(r'\toprule')
        ncols = len(types + [t for t in types if t == 1]) # Single err case
        self.header.append(r'{} & \multicolumn{'+str(ncols-1)+'}{c}{'+title+r'} \\')
        self.header.append(r'\cmidrule(l){2-'+str(ncols)+'}')

        # Labels for units
        lablst = []
        for t, lab in zip(types, labs):
            if t == 1:
                lablst.append(r'\multicolumn{2}{c}{' + lab + r'} & ')
            else:
                lablst.append('{} & '.format(lab))
        # Remove trailing ampersand, add newline
        self.header.append((''.join(lablst))[:-3] + r'\\')
        self.header.append(r'\midrule')

        # Build row spec
        rlst = []
        for t in types:
            if t == 0:
                rlst.append('{0.3g} & ')  # Default for a number
            elif t == 1:
                rlst.append('${}$ & ${}$ & ')  # Single error
            elif t == 2:
                rlst.append('${{{}}}^{{{}}}_{{{}}}$ & ')  # Double error
                # Need braces around data value, in exponential case
            else:
                rlst.append(t + ' & ')  # Supply your own damn string spec
        # Remove trailing ampersand, add newline
        self.rspec = (''.join(rlst))[:-3] + r'\\'

        # Rows to be filled with rspec-formatted strings
        self.rows = []

        self.footer = [r'\bottomrule', r'\end{tabular}']


    def __str__(self):
        return '\n'.join(self.get_table())

    def get_table(self):
        return self.header + self.rows + self.footer

    def add_row(self, *args):
        """Add row to table, using rspec.
        If giving two errors, positive err comes first..."""
        args = list(args)  # Must expand tuple
        ind = 0
        for t in self.types:
            if t == 1:
                args[ind:ind+2] = self.fmt_numerr(*args[ind:ind+2])
                ind += 2
            elif t == 2:
                args[ind:ind+3] = self.fmt_num_2err(*args[ind:ind+3])
                ind += 3
            else:
                ind += 1
        self.rows.append(self.rspec.format(*args))

    # TODO this code is very gross
    # Sept. 10 -- forget the precision business.  Let the user deal with
    # proper rounding, report more than needed in floating pt form

    def fmt_num_2err(self, num, errpos, errneg):
        nstr = '{:0.3g}'.format(num)
        estr_pos = '+{:0.2g}'.format(errpos)
        estr_neg = '-{:0.2g}'.format(errneg)

        return self.strfl(nstr), self.strfl(estr_pos), self.strfl(estr_neg)
    
    def fmt_numerr(self, num, err):
        nstr = '{:0.3g}'.format(num)
        estr = '{:0.2g}'.format(err)
        return self.strfl(nstr), self.strfl(estr)

    def strfl(self, x):  # Yes, really.  See http://stackoverflow.com/q/3410976
        return str(float(x))

#    def fmt_num_2err(self, num, errpos, errneg):
#        """Format number w/error nicely.
#        It's imperfect and doesn't agree w/ rounding of Ressler,
#        but it's simpler this way...
#        DEAR GOD THIS CODE IS DISGUSTING
#        """
#        num_exp = int(np.floor(np.log10(num))) if num != 0. else 0
#        err_min = abs(min(errpos, errneg))
#        err_exp = int(np.floor(np.log10(err_min))) if err_min != 0. else 0
#
#        prec = num_exp - err_exp + 2 # Print to precision + 1, allowed by error
#        if prec <= 0:
#            prec = 1  # Enforce minimum precision
#
#        # Use string formatting to handle rounding
#        numstr = ('{:0.'+str(prec)+'g}').format(num)
#        errstr_pos = ('+{:0.2g}').format(errpos)
#        errstr_neg = ('-{:0.2g}').format(errneg)
#
#        numstr = self.expstr2latex(numstr)
#        errstr_pos = self.expstr2latex(errstr_pos)
#        errstr_neg = self.expstr2latex(errstr_neg)
#        return numstr, errstr_pos, errstr_neg
#
#    def fmt_numerr(self, num, err):
#        """Format number w/error nicely.
#        It's imperfect and doesn't agree w/ rounding of Ressler,
#        but it's simpler this way...
#        """
#        num_exp = int(np.floor(np.log10(num))) if num != 0. else 0
#        err_exp = int(np.floor(np.log10(err))) if err != 0. else 0
#        prec = num_exp - err_exp + 2 # Print to precision + 1, allowed by error
#        if prec <= 0:
#            prec = 1  # Enforce minimum precision
#
#        # Use string formatting to handle rounding
#        numstr = ('{:0.'+str(prec)+'g}').format(num)
#        errstr = ('{:0.2g}').format(err)
#
#        numstr = self.expstr2latex(numstr)
#        errstr = self.expstr2latex(errstr)
#        return numstr, errstr

    def expstr2latex(self, numstr):
        """Convert a Python string-formatted number to nicer LaTeX"""
        # Capturing groups for significand, exp sign if negative, exp power
        exp_regexp = r'([0-9\+\-\.]+)' + 'e' + r'([\-]*)[\+0]*([1-9]+[0-9]*)'
        result = re.search(exp_regexp, numstr)
        if result:
            return (r'{} \times '.format(result.group(1)) +
                    '10^{' + '{}{}'.format(*result.group(2,3)) + '}')
        return numstr


if __name__=='__main__':
    pass

