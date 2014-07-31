Log of tables and generating parameters/grids
=============================================

Tables of FWHM values over parameter space of (mu, eta2, B0), generated from
Sean Ressler's full model code.  Tables and associated log files are READ-ONLY
(currently just running chmod a-w [filename] manually), to avoid overwriting
tables for now.

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

`sn1006_grid-1-100-20_2014-07-28_mu=0.25_partest.pkl`
---------------------------------------------
And, three additional variants of this (mu=0.66, mu=0.75, mu=1.25)
Checking Brian's ad hoc parallelization idea.  With 4 cores, all four of these
should run in approx 1 hour!  Start them all at abt the same time.

    mu_vals = [mu]
    eta2_vals = np.logspace(-2, 2, 50, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data/1.51, axis=0)
    data_max = np.amax(data*1.51, axis=0)

    fname = 'sn1006_grid-1-100-20_2014-07-28_mu-{:0.2f}_partest.pkl'.format(mu)

    # Buried inside code (kind of, now included in log file)
    rminarc = max(fwhms_max) * f_minarc
    f_minarc = 1.2
    f_B0_init = 1.1
    f_B0_step = 0.15

All started at ~15:01.  Yup, looks like this is taking full advantage of
all four cores!  Yay!
Finished: 19:16, 19:16, 18:12, 18:26.  Interesting.

Compare -- 6 mu values, w/ otherwise same grid.  Start 22:51, finish 17:30 the
next day (total = 17:30 + 1:09 = 18:09, means ~3 hrs / mu value).

Here we finished at 19:08 on average --> 4:07 hrs / 4 mu values, or about 1 hr/
mu value now.  So we cut the time expenditure by 1/3rd!


`tycho_grid-6-100-20_...`


Changed  in `maketab_gridB0`, dx_max (i.e. the B0 step) to 15 percent of the
initial B0 value.  Yes, it's another twiddle-able...

Gridding for this on hold while I address rminarc problem.
