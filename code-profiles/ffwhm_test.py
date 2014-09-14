"""
Really simple tests for FWHM calculation and stretching
From command line, run nosetests

Originally run inside IPython
Using [Taavi's IPython nose](https://github.com/taavi/ipython_nose).
Some useful links:
* [Nose intro](http://pythontesting.net/framework/nose/nose-introduction/)
* [unittest.TestCase API](https://docs.python.org/2/library/
unittest.html#unittest.TestCase) (same stuff in nose.tools)

Aaron Tran
July 2014; split to .py file Sept 2014
"""

import numpy as np
from nose.tools import assert_almost_equal, assert_equal, assert_true

import ffwhm

def gaussian(x, x0, sigma):
    """Not normalized, just for testing.
    mathworld.wolfram.com/GaussianFunction.html
    Expected FWHM is 2*sqrt(2 ln(2))*sigma
    endpoints at x0 +/- sqrt(2 ln(2))*sigma (~ x0 +/- 1.18*sigma)
    """
    return np.exp( -(x-x0)**2 / (2*sigma**2) )


def test_FWHM_gaussian():
    """Check FWHMs of various gaussians"""

    x0_vals = [-2., 0., 1., 1e4]
    sigma_vals = [0.02, 1, 1.5, 3.]

    for x0, sigma in zip(x0_vals, sigma_vals):
        fwhm_ans = 2*np.sqrt(2*np.log(2))*sigma

        # x_max not specified
        fwhm, fwhm_a, fwhm_b = ffwhm.get_fwhm(x0-1.2*sigma, x0+1.2*sigma, gaussian, [x0, sigma])

        assert_almost_equal(fwhm, fwhm_ans)
        assert_almost_equal(fwhm_a, x0 - fwhm_ans/2.)
        assert_almost_equal(fwhm_b, x0 + fwhm_ans/2.)

        # x_max not specified, with crappy search bounds
        fwhm, fwhm_a, fwhm_b = ffwhm.get_fwhm(x0-100*sigma, x0+1.1776*sigma, gaussian, [x0, sigma])

        assert_almost_equal(fwhm, fwhm_ans)
        assert_almost_equal(fwhm_a, x0 - fwhm_ans/2.)
        assert_almost_equal(fwhm_b, x0 + fwhm_ans/2.)

        # x_max specified
        fwhm, fwhm_a, fwhm_b = ffwhm.get_fwhm(x0-2*sigma, x0+2*sigma, gaussian, [x0, sigma], x_max=x0)

        assert_almost_equal(fwhm, fwhm_ans)
        assert_almost_equal(fwhm_a, x0 - fwhm_ans/2.)
        assert_almost_equal(fwhm_b, x0 + fwhm_ans/2.)

        # Does limit checking work?
        fwhm, fwhm_a, fwhm_b = ffwhm.get_fwhm(-1*sigma, 1*sigma, gaussian, [x0, sigma])
        assert_true(np.isnan(fwhm))


def test_stretch_function():
    """Simple checks for stretch function (not comprehensive)"""

    # No stretch occurs for xi = 0
    assert_equal(ffwhm.stretch_x(-15, 0, 10), -15)
    assert_equal(ffwhm.stretch_x(7, 0, 10), 7)

    xi_vals = [1, 2, 3]
    x0_vals = [1, 10, 100, -10]
    for xi, x0 in zip(xi_vals, x0_vals):
        assert_equal(ffwhm.stretch_x(0, xi, x0), 0)  # Zero invariant
        assert_equal(ffwhm.stretch_x(x0, xi, x0), x0)  # x0 invariant


def test_inv_stretch_function():
    """Simple checks for inverse stretch function (not comprehensive)"""
    xi = 1.
    x0 = 25
    xc = 55

    for x in np.linspace(0, 55, 13):  # Check over range [0, 55]
        assert_almost_equal(x, ffwhm.inv_stretch_x(ffwhm.stretch_x(x,xi,x0,xc), xi,x0,xc))
        assert_almost_equal(x, ffwhm.stretch_x(ffwhm.inv_stretch_x(x,xi,x0,xc), xi,x0,xc))


def test_FWHM_stretched_gaussian():
    """Check FWHMs of stretched gaussians"""

    # NB this is sensitive to values of xi, x0_stretch, etc.
    # MUST choose sensible values or things will blow up.
    # Note that stretching is generally asymmetric

    # First, get unstretched Gaussian fwhm and fwhm bounds
    x0_gauss = 20
    sigma = 5
    fwhm, a, b = ffwhm.get_fwhm(x0_gauss-2*sigma, x0_gauss+2*sigma, gaussian, [x0_gauss, sigma])
    print fwhm, a, b

    # Stretch parameters
    xi = -0.4  # Negative xi widens the Gaussian
    x0_stretch = x0_gauss  # Must always be the case

    # Compute new FWHM by stretching old FWHM interval
    a_str = ffwhm.inv_stretch_x(a, xi, x0_stretch)
    b_str = ffwhm.inv_stretch_x(b, xi, x0_stretch)
    fwhm_str = b_str - a_str
    print fwhm_str, a_str, b_str

    # Compare against stretched function
    gaussian_str = ffwhm.stretch_func(gaussian, xi, x0_stretch)

    assert_almost_equal(0.5, gaussian_str(a_str, x0_gauss, sigma))
    assert_almost_equal(0.5, gaussian_str(b_str, x0_gauss, sigma))

if __name__ == '__main__':
    pass
