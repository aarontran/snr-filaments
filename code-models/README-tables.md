Log of tables and generating parameters/grids
=============================================

Current (ad hoc) naming convention is:
snr name, grid resolution (3-tuple) of mu, eta2, B0 (in that order),
date (yyyy-mm-dd).

`sn1006_grid-3-40-15_2014-07-25.pkl`
------------------------------------

A first attempt, after one failure in the morning.
The data are good, but I did not sample widely enough in FWHM space

    mu_vals = np.array([1/3, 1/2, 1])  # Kolmogorov, Kraichnan, Bohm
    eta2_vals = np.logspace(-2, 2, 40, base=10)
    n_B0 = 15
    
    data_min = np.amin(data, axis=0)
    data_max = np.amax(data, axis=0)
    
    # Buried inside code... you should let this be a free setting, perhaps...
    rminarc = max(fwhms_max)*f_minarc
    f_minarc = 1.2

I'm not certain of the value of `f_minarc` here, may have been 1.1 or 1.2.
But, I think it was 1.2.


`sn1006_grid-6-100-20_2014-07-25.pkl`
-------------------------------------
(along with associated .log, .errlog)

6 x 100 x 30 x (3 s) = 15 hrs, nice and solid.
To avoid redundancy in eta2 values, I take 50 logarithmically spaced points
and 50 linearly spaced points, following Brian's suggestion.

This is important when best fit eta2 is not ~0.

    mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]
    eta2_vals = np.logspace(-2, 2, 50, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20
    
    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data/1.51, axis=0)
    data_max = np.amax(data*1.51, axis=0)
    
    # Buried inside code
    rminarc = max(fwhms_max)*f_minarc)
    f_minarc = 1.2

The code appears to go bonkers for very specific values of rminarc.
Hence the use of scaling factor 1.51 to determine `data_min` and `data_max`.
See my debugging notes on this matter...