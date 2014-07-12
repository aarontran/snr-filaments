"""
Small module with methods for matplotlib
Goal is to make plotting easier, and ensure reproducible plots
(and, cut down on amt of boilerplate code)

fplot is intended if you are using the interactive like interface,
rather than the object-oriented interface.  A convenience method only.

This is supposed to be small enough to be okay for global import
from fplot import *

Aaron Tran
6 June 2014
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

def main():
    show_mplrc_settings()


def show_mplrc_settings():
    """Display information about matplotlibrc file"""
    print 'Using %s' % mpl.matplotlib_fname()
    r = mpl.rcParams

    ff = r['font.family'][0]
    print 'Font sizes for axes: %g; (x,y) ticks: (%g, %g): legend %g' % \
          (r['axes.labelsize'], r['xtick.labelsize'],
           r['ytick.labelsize'], r['legend.fontsize'])
    print 'Font family %s uses face %s' % (ff, r['font.'+ff])

    print 'Figure size: %s, dpi: %g' % (r['figure.figsize'], r['figure.dpi'])


def fplot(xlab, ylab, ax=None, axargs=None, scales=None):
    """Convenience method to set x/y labels, axis limits in one cmd
    xlab: xlabel string
    ylab: ylabel string
    axargs: limits [xmin, xmax, ymin, ymax] or other string arg (e.g., 'tight')
    scales: [xscale, yscale] (see plt.xscale(...))
    """
    if ax is None:
        ax = plt.gca()
    if scales is not None:  # Could allow to accept only one arg
        ax.set_xscale(scales[0])
        ax.set_yscale(scales[1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if axargs is not None:
        ax.axis(axargs)
    # plt.ticklabel_format(useOffset=False, axis='x')
    # plt.ticklabel_format(useOffset=False, axis='y')


if __name__ == '__main__':
    main()

