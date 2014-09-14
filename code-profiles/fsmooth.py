"""
Functions to smooth profiles and identify extrema for fit domains.
The first minimum downstream of the thin rim bounds the fit domain.

Aaron Tran
July 2014 (split to .py file, Sept 2014)
"""

import numpy as np
import scipy as sp
from scipy import optimize
from scipy import signal


def ind_first_min(y, **kwargs):
    """Index of first minimum to left of global non-edge max in smoothed data
    If no such minimum is found, return 0 (i.e., don't cut the fit domain).
    If no max (peak) is found, return 0.
    **kwargs (window_len, window) passed to std_smooth
    """

    # List of indices of relative minima in smoothed-y
    y_smooth = std_smooth(y, **kwargs)
    inds_relmin = sp.signal.argrelmin(y_smooth)[0]  # Unpack the tuple...

    if len(inds_relmin) is 0:
        return 0

    # Index of global (non-edge) maximum in y
    try:
        ind_max = get_ind_max(y, **kwargs)
    except ValueError:
        return 0

    # Find the largest value in inds_relmin, that's smaller than ind_max
    max_left_min = 0
    for ind_min in inds_relmin:
        if ind_min > max_left_min and ind_min < ind_max:
            max_left_min = ind_min

    return max_left_min

def get_ind_max(y, **kwargs):
    """Smooths data and returns index of largest relative maximum in smoothed data
    **kwargs (window_len, window) passed to std_smooth

    If no relative maximum found,
    """
    ysmooth = std_smooth(y, **kwargs)
    # Get indices of relative maxima
    inds_relmax = sp.signal.argrelmax(ysmooth)[0]  # Unpack the tuple...
    # Get index of the largest relative maximum
    if len(inds_relmax) == 0:
        raise ValueError("No relative maximum found in y for given smoothing")

    idx = np.argmax(ysmooth[inds_relmax])  # idx indexes inds_relmax, now
    return inds_relmax[idx]

def std_smooth(x, window_len=21, window='hanning'):
    """Modified from scipy cookbook (http://wiki.scipy.org/Cookbook/SignalSmooth)
    Add reflected copies of signals at each end to address convolution edge effects

    Default window length of 21 pts, empirically works well for Tycho's SNR.
    No problems with hanning, haven't tested other windows.  Choice shouldn't matter much.

    Input:
        x: array to be smoothed
        window_len, window: window length/type.  Best if window_len is odd,
                            but works for even window length (haven't checked
                            exactly how np.convolve handles even convolution,
                            so I just used mode='same' to get symmetric output.
    Output:
        y: smoothed array of same size as x (padding for convolution removed)
    """
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    pad_len = window_len - 1

    if window == 'flat': #  moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.{}(window_len)'.format(window))

    y = np.convolve(w/w.sum(), s, mode='same')
    return y[pad_len:-pad_len]


if __name__ == '__main__':
    pass

