"""
Currently, single class definition of a LatexTable
for formatting numpy data with errors into a nice LaTeX'd table.
Designed for LaTeX package 'booktabs'

Aaron Tran
August 2014


This duplicates a LOT of functionality from `matrix2latex.py`.
I'd like to extend matrix2latex by just adding support for data w/errors;
it handles a lot of other functionality nicely already.
E.g., see the implementation of niceFloat
* use errors to compute # of sigfigs and supply that as an argument
* use \e in LaTeX to change exponential notation formatting, rather than
  manually editing floats


"""

from __future__ import division
import numpy as np
import re

#from matrix2latex import matrix2latex as m2ltx

# Example matrix2latex implementation:
# algnmt = '@{}rr@{ $\pm$ }lr@{ $\pm$ }lr@{}'
# hr = [[''] + 5 * [title],
#       [r'$\mu$ (-)'] + 2 * [r'$\eta_2$ (-)] + 2 * [r'$B_0$ ($\mu$G)']
#        + [r'$\chi^2$']]
#
# fmtcol = 'STUFF'  # Here parse sigfigs, depending on error
# # PROBLEM... fmtcol will vary depending on the row.
# # So this is where I probably need to exted matrix2latex
# m2ltx(data, headerRow = hr, formatColumn = fmtcol,
#       alignment = algnmt)


class LatexTable(object):
    """Worth generalizing, someday.
    1. customize number of columns, +/- formatting, number formatting
    2. customize concatenating tables together... (horizontally, that is)
    3. customize header/footer/layout control

    Setup for booktabs table output (what about aastex deluxetable?)

    Options for columns:
    1. number w/ single +/- error
        tablespec: + 'r@{ $\pm$ }l'
        row fmt: + ' & ${}$ & ${}$'
        units header: + r'\multicolumn{2}{c} ...'
    2. number w/ distinct +/- errors
        tablespec: + 'r'
        row fmt: + ' & ${}^{}_{}$'
        precision of stated value set by smallest error value.
    3. something else, e.g. single # (give the fmt specifier directly)
        tablespec: + 'r'
        row fmt: + ' & STUFF' (you decide precision)

    add_row then accepts an arbitrary number of columns
    """

    # Initializes: header, labels, rows, footer; rspec, types
    def __init__(self, labs, types, title):
        """Types specify column types; 0,1,2 = no/single/double +/- errors"""

        self.types = types

        # Build table/column spec, treating t==1 (single err) case specially
        hlst = [r'\begin{tabular}{@{}']
        for t in types:
            if t == 1:
                hlst.append(r'r@{ $\pm$ }l')
            elif t == 2:
                hlst.append('l')
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
    
    #def add_row(self, mu, e2, e2eps, B0, B0eps, chisqr):
    #    """Add row to table, w/assumed format"""
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
        #self.rows.append(fmtstr.format(mu, e2_str, e2eps_str, B0_str,
        #                               B0eps_str, chisqr))

    def fmt_num_2err(self, num, errpos, errneg):
        """Format number w/error nicely.
        It's imperfect and doesn't agree w/ rounding of Ressler,
        but it's simpler this way...
        DEAR GOD THIS CODE IS DISGUSTING
        """
        num_exp = int(np.floor(np.log10(num))) if num != 0. else 0
        err_min = abs(min(errpos, errneg))
        err_exp = int(np.floor(np.log10(err_min))) if err_min != 0. else 0

        prec = num_exp - err_exp + 2 # Print to precision + 1, allowed by error
        if prec <= 0:
            prec = 1  # Enforce minimum precision

        # Use string formatting to handle rounding
        numstr = ('{:0.'+str(prec)+'g}').format(num)
        errstr_pos = ('+{:0.2g}').format(errpos)
        errstr_neg = ('-{:0.2g}').format(errneg)

        numstr = self.expstr2latex(numstr)
        errstr_pos = self.expstr2latex(errstr_pos)
        errstr_neg = self.expstr2latex(errstr_neg)
        return numstr, errstr_pos, errstr_neg

    def fmt_numerr(self, num, err):
        """Format number w/error nicely.
        It's imperfect and doesn't agree w/ rounding of Ressler,
        but it's simpler this way...
        """
        num_exp = int(np.floor(np.log10(num))) if num != 0. else 0
        err_exp = int(np.floor(np.log10(err))) if err != 0. else 0
        prec = num_exp - err_exp + 2 # Print to precision + 1, allowed by error
        if prec <= 0:
            prec = 1  # Enforce minimum precision

        # Use string formatting to handle rounding
        numstr = ('{:0.'+str(prec)+'g}').format(num)
        errstr = ('{:0.2g}').format(err)

        numstr = self.expstr2latex(numstr)
        errstr = self.expstr2latex(errstr)
        return numstr, errstr

    def expstr2latex(self, numstr):
        """Convert a Python string-formatted number to nicer LaTeX"""
        # Capturing groups for significand, exp sign if negative, exp power
        exp_regexp = r'([0-9\+\-\.]+)' + 'e' + r'([\-]*)[\+0]*([1-9]+[0-9]*)'
        result = re.search(exp_regexp, numstr)
        if result:
            return (r'{} \times '.format(result.group(1)) +
                    '10^{' + '{}{}'.format(*result.group(2,3)) + '}')
        return numstr


if __name__ == '__main__':
    pass

