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


def fplot(xlab, ylab, axargs='tight', scales=('linear','linear')):
    """Convenience method to set x/y labels, axis limits in one cmd
    xlab: xlabel string
    ylab: ylabel string
    axargs: limits [xmin, xmax, ymin, ymax] or other string arg
    scales: [xscale, yscale] (see plt.xscale(...))
    """
    plt.xscale(scales[0])
    plt.yscale(scales[1])
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.axis(axargs)
    # plt.ticklabel_format(useOffset=False, axis='x')
    # plt.ticklabel_format(useOffset=False, axis='y')


if __name__ == '__main__':
    main()

