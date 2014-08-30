Research notes
==============

Aaron Tran
Fall 2014

(continued from rsch-notes.md)

Timeline
========

Weeks 11-13: late August
Weeks 14-18: September (AAS abstract deadline is week 18)
Weeks 18-22: October
Weeks 23-26: November (Thanksgiving, week 26)
Weeks 27-30: December (Christmas, week 30)


Table of contents
=================

* Week 10 - full model grid-best-"fits" (test fits, error annealing w/ SN 1006)
            port full model code to Python, optimize for speed
* Week 11 - one week break
* Week 12 - code for LaTeX tables; more precise manual error calculations.
            Debug/test new full model
* Week 13 - refactor model exec/disp code.  Debug error calculations
            extensively.

(week 10 included for continuity)


(Week 11) Monday 2014 August 11 -- Friday 2014 August 15
========================================================

On vacation


(Week 12) Monday 2014 August 18
=====================

Summary
-------
Tables...

Set-up for mac work
-------------------
For now, until I get a computer/workstation at GSFC again, work off of my mac's
local git repo.

Trying to work with `fullmodel.so` (f2py compiled model code)
* Problems: no fortran compiler found (f2py output).  Couldn't locate
libgfortran.3.dylib in `/opt/local/lib/gcc/` or something...
* Solution: sudo port select --set gcc mp-gcc48

As of TODAY? my access to Berkeley's library proxy is dead.
Use UMD Research Port now...


Table generation, errors
------------------------
Major issue: method of confidence interval contours breaks down at large
chisquared... but I don't know when/how, or which references to seek.
I thought the `lmfit` package had some info on this but I can't find it.
Perhaps see [this](http://stats.stackexchange.com/q/76444).
And, [this](http://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm).
Anyways...


Tuesday 2014 August 19
======================

Summary
-------
* More robust FWHM computation in full model code (Python port) addresses
  sharp drop-off at r=1 (box length error; fall-off at rim edge is weird)
* Fail more gracefully if FWHM is not found (extreme values of B0, eta2)
* Fix FWHM calculation bugs at intensity grid edge (for very narrow rims)


Result: rminarc must be sized to fit rims still, AND it should be as small as
possible to reduce resolution error (especially at high energies).  Quantifying
that comes soon.

Remark: wow, I have been working on the "full model code" for some 4 weeks
now... minus a bit of time for the poster presentation.

Learn to program more effectively and avoid bugs in the first place.
Because time is expensive...  Think preventative coding.


Remove box length errors at r=1
-------------------------------
A bug appears on a full model call with:

    B0 = 821.815 muG; eta2 = 9.796; mu = 0.000;
    rminarc = [ 18.5    15.    10.77  10.77  10.77]

We get box length errors (max r-value in rmesh is not close enough to edge, to
capture rim drop off).  Output FWHMs are quite small:

    [ 1.9425    1.725     1.373175  1.1847    1.02315]

I fixed this bug; code now searches for FWHM crossings on right.


Bug in SN 1006 test (rim too smeared out)
-----------------------------------------

Code fails on function call with (using SN 1006 parameters):

    B0 = 67.933 muG; eta2 = 1187.569; mu = 0.000; rminarc = [60.  60.  60.]

Code fails while finding the intensity maximum; max is found at smallest
rmesh value (r = rmin).  It seems that B0 was too small, eta2 was too big, so
the rim was smeared out horribly.

Solution: print warning and return error value (fwhm = shock radius).
This allows fitting to proceed if we hit extreme parameter values.


rminarc-dependent FWHMs at high energy
--------------------------------------

Result: bugs in intensity-finding method.

Weird issue -- at large energy (10 keV), rminarc strongly affects FWHM values.
All subsequent debugging uses the following parameters:

    B0 = 300e-6, eta2 = 1000, mu = 0; Tycho's SNR with vs = 4.52e8

At 10 keV, changing rminarc from 10 to 200, changes output FWHM by ~0.25 arcsec
(for true FWHM ~6.95 arcsec).  This 4% error is quite large!
(note that for 0.7 keV the error remains small, looks okay)

### Debugging log

Changing rminarc changes the peak intensity value a lot!

    peak intensity = 6.31849924787e-135 at r = 0.994645003207 (rminarc = 10.)
    peak intensity = 7.35531230369e-135 at r = 0.994283560673 (rminarc = 200.)

Hypothesis: is scipy minimization to find peak intensity, not consistent?
I was wondering if it would depend on the epsilon tolerance.
At rminarc = 10, 200, I computed intensities at the two different peak
locations (r = 0.9946450, 0.9942835).
Result: computed intensity changes a lot, whether on/off peak, depending on
rminarc.  Scipy finds the peak correctly in each case.

Changing rminarc, then, is changing overall intensity values.  I'm befuddled.
Why does rminarc change intensity calculation?  Looking in code:
* `rmin_nu` sets values of e- distr to interpolate over (rmin to 1)
  (sets resolution of disttab)
* `rmin_nu` also sets values of emissivity to interpolate over (rmin to 1)
  (sets resolution of emistab)
These two bullet points are the same, as emissivity just integrates e-
distribution over electron energy.

If I change iradmax to improve e- distr grid, does FWHM w/rminarc=200 get
closer to FWHM w/rminarc=10? (assuming that smaller rminarc means better
resolution for e- distr interpolation)
CONFIRMED: iradmax = 100 vs. 400 makes a HUGE difference.
With iradmax=400, the FWHMs agree to 0.007 arcsec (vs. 0.25 arcsec before).
Absolute FWHM value changed by ~0.0003 for rminarc = 10 arcsec.

(sanity check) changing irmax from 100 to 400 (holding iradmax=100, ixmax=500)
made no difference, as expected.

### Are rminarc, iradmax interchangeable in setting resolution?

Simple test: does doubling iradmax have the same effect as halving rminarc?
Same parameters as before (irmax=100, ixmax=500; 10 keV energy w/ Tycho's SNR)

    iradmax = 100, rminarc = 20: FWHM = 6.9568692 (baseline)
    iradmax = 100, rminarc = 10: FWHM = 6.95696361 (change = 0.00009441)
    iradmax = 200, rminarc = 20: FWHM = 6.95697527 (change = 0.00010607)

So, doubling iradmax / halving rminarc effects a change of ~1e-4, but each
method differs by ~1e-5.  Looks promising.

How about the effects of a larger/sharper change?  Do the methods differ by
about 10% normally, or do they differ by an absolute amount ~1e-5 arcsec?

    iradmax = 100, rminarc = 100: FWHM = 7.15124535
    iradmax = 200, rminarc = 200: FWHM = 6.97627089

Goddamnit (difference: ~0.18).  Maybe that shouldn't be a surprise.
Repeat this for 0.7 keV and 4.5 keV:

    iradmax = 100, rminarc = 100: FWHM(0.7, 4.5 keV) = 26.21491971, 10.35922423
    iradmax = 200, rminarc = 200: FWHM(0.7, 4.5 keV) = 26.21539187, 10.3612368

Differences are ~0.0005, 0.002.  Wow.  The discrepancy increases very sharply
with energy.  For comparison, check against "good" values.

    iradmax = 500, rminarc = [30, 12]: FWHMs = 26.22028134, 10.37183778
    iradmax = 900, rminarc = [30, 12]: FWHMs = 26.22013231, 10.37182883

Cutting down rminarc / increasing iradmax, gives changes of ~0.005, ~0.01.


### Apparently not.  What's going on?

Questions:
1. why are these effects not interchangeable?
2. why does the resolution error worsen with larger energy? (and so sharply, it
   almost looks exponential)

Note that: d(rad) = (1 - rmin)/(iradmax - 1) = rminarc/rsarc /(iradmax - 1)
To be fair, the comparison should "strictly" be (at 10 keV):

    iradmax = 100, rminarc = 100: FWHM = 7.15124535
    iradmax = 199, rminarc = 200: FWHM = 7.47571061
    iradmax = 200, rminarc = 200: FWHM = 6.97627089

Wait, what?!  This deeply bothers me not because it affects our results much,
but because this means I don't truly understand the code / model...

I walk through code and try to find where values diverge, with iradmax/rminarc
set to 100/100 and 199/200, which I expect to give same results.
CHECK: e- and emissivity distributions (disttab, emistab) appear to agree,
from casual inspection, plots of small subsets.

What happens if I update irhomax to match our scaling of rminarc?
From visual inspection of a tiny subset, with
iradmax/rminarc/irhomax = 100,100,10000 and 199,200,19999,
emistab is identical.  Let's calculate FWHMs

    iradmax = 100, rminarc = 100, irhomax = 10000: FWHM = 7.1512
    iradmax = 199, rminarc = 200, irhomax = 10000: FWHM = 7.47571
    iradmax = 199, rminarc = 200, irhomax = 19999: FWHM = 6.9784  (?!?!?!?!)

I know the FWHM is wrong, but I expect it to be _consistently_ wrong, damnit.

Since disttab, emistab are the same, I expect intensity computed at a given
r-value to give the same result.  Checking this: AFFIRMED.
The code with iradmax/rminarc = 100/100 is NOT finding the local max.

CHECK: reenable `xatol` for SciPy's `minimize_scalar`.  Nothing happens.

### Realization -- I'm hitting the r-grid (irmax) edge on the RIGHT side.

I think I realized the problem. At the grid's right edge (`rmesh[-1]`),
I seek a local maximum between rmesh[-1], 1.
But, the maximum may lie between rmesh[-2], rmesh[-1].

Example:

    rminarc/iradmax/irhomax = 100/100/10000

    Searched between [0.995833, 0.997917] (intensities: [6.415, 5.600] x1e-135)
    Found max at r = 0.995833351075 (left edge was 0.995833333333)

    rminarc/iradmax/irhomax = 200/199/19999
    (didn't change irmax, so coarser intensity grid)

    Searched between [0.991667, 0.997917] (intensities: [6.235, 5.600] x1e-135)
    Found max at r = 0.99419 < 0.995833 (left bound for rminarc=100 case!)

And, as expected, the rminarc=200 case got a more accurate FWHM.

PATCH THE BUG.  That seems to fix things.  Let's run some tests

    rminarc iradmax irhomax FWHM

    100     100     10000   6.97837734
    200     199     19999   6.97837734  YAY!
    200     200     20000   6.97628278  Looks okay

    200     199     10000   7.47571707  WHAT?!
    200     199     20000   7.47566551  ?!?!?!

    ... (drastic increase in resolution)

    20      100     10000   6.9568692
    20      199     10000   6.95696459
    10      100     10000   6.95696359 (good, irhomax doesn't matter)
    10      500     20000   6.95710649

ONE MORE PATCH -- floating point error, when searching for a crossing, take a
step first them compare intensity values.  In the weird cases (fwhm=7.48) code
didn't search for a crossing, hence the error.

    200     199     10000   6.97832483  Now this looks okay too
    200     199     20000   6.97837896

Please note that these bugs appear ONLY if irmax is too small to
resolve intensity peak on the right (i.e., peak is too narrow).

And, this should explain why the error was so large at high energy --
at high energy we are hitting the edge of the irmax grid, and hence running
into these bugs.

### Parting remarks

In practice, I wouldn't have caught these errors as they crop up rarely
(at lower energy, larger irmax, smaller rminarc).
But it does have the potential to mess up results.  Very subtle bugs...

I'm afraid, what other bugs are there that I don't know about.
As Michael reminds me: write unit tests, then write your code...

This also underscores the need for a CAREFUL validation of our model over
parameter space.  We can't blithely leave iradmax, ixmax, rminarc at arbitrary
default settings.


Wednesday 2014 August 20
========================

Summary
-------
* Removed old B0 interpolation code
* Add statement to recompile f2py code in Python full model port
* Code tidying, docstring updates, etc
* Error annealing now 1. uses rootfinding for errors, 2. scans outside grid


Note: code folding in Vim is really handy.  Makes it a lot easier to see
1. method hierarchy, inputs/outputs/1-liner docstrings
2. unfolded code indicates work in progress (updated .vimrc to save code
   folding state between sessions)


B0 interpolation code
---------------------
Reviewed and removed from Tycho ipynb.  Interpolation still seems uncertain and
inconsistent; it doesn't consistently find better chisqr values (using the 
optimized, more accurate/faster Python port of full model code).
Note that B0 interp code was removed in commit msg.

This may be because the table FWHMs are from old full model code (Sean's), not
the new Python port (which should vary continuously and be more accurate).

If you want to retry this, generate new tables and yank code from old commit.
If interpolation gives consistently better B0 values it may be worth it.
For now, it complicates things, and may be covered by fitting (free eta2, B0)
from best grid value anyways.


Testing code on SN 1006, once more
----------------------------------
Now that (more) bugs are addressed, code seems to run smoothly.
I run the code for 20--40 minutes on SN 1006, Filament 1 w/ 3 bands.

(full model fits try several fits w/ multiple epsfcn values, starting from
grid best fit values)

----------------------------------------------------------------

### Results (full model fits for params, and error annealing):

    mu      eta2                    B0                      chisqr
    0.00    27.998 (+72.002/-26.672)187.92 (+55.86/-79.98)  0.1149
    0.33    5.510 (+94.490/-4.490)  126.67 (+104.91/-26.29) 0.0173
    0.50    4.286 (+95.714/-3.375)  117.82 (+105.79/-20.14) 0.0399
    1.00    2.589 (+5.778/-1.679)   102.49 (+18.51/-8.67)   0.1415
    1.50    2.458 (+3.869/-1.547)   97.08 (+10.58/-5.85)    0.2221
    2.00    2.640 (+3.074/-1.620)   94.18 (+6.44/-4.34)     0.3021

For mu=0.00, 0.33, 0.50, eta2's upper limit = 100.0 (due to eta2 grid).
Next, compare to best _grid_ values (I believe errors are not annealed)

### Results (best grid values [full model])

    mu      eta2                    B0                      chisqr  change
    0.00    10.000 (+58.665/-8.367) 151.46 (+75.70/-40.78)  0.1748  +0.06
    0.33    5.510 (+77.354/-4.286)  126.67 (+97.65/-24.84)  0.0173  n/a?
    0.50    4.286 (+78.579/-3.187)  117.82 (+92.30/-18.80)  0.0399  n/a?
    1.00    2.857 (+5.510/-1.947)   102.98 (+19.00/-9.94)   0.1973  +0.06
    1.50    1.931 (+4.396/-0.910)   94.71 (+13.84/-2.84)    0.3516  +0.13
    2.00    2.330 (+3.180/-1.004)   93.83 (+5.37/-1.92)     0.3750  +0.07

Why did fits not improve from grid values at mu=0.33, 0.50?
It seems unlikely that the "best grid values" would be a local minimum in
chisqr space exactly.

Hypothesis: because I changed the full model code, the tabulated FWHMs may be
erroneously "good".  So, when running full model fit, we can't find any better
FWHMs than the pre-cached ones (even at the exact same data points).

Need to rerun tables _and_ recompute FWHMs at identified best grid values.
Use tabulated FWHMs to find best grid (eta2, B0), but always recompute FWHMs at
this "best" value.

### Time analysis

For mu = 0.00 case:
* 121 function calls to find best fit, started from best grid values
* 42 function calls to anneal errors on left
* 28 function calls to anneal errors on right (but, I hit grid edge)
* Total: 191 function calls

200 calls x (6 mu values) x (3 seconds) ~ 1 hr, for 1 region/filament.
(may need more function calls for better error estimates)

So ~5+ hrs for SN 1006; ~13-18+ hrs for Tycho's SNR.

### Conclusions, next steps

1. error annealing should optimize to find chisqr threshold and scan outside of
   grid values. Currently, we anneal each grid point until we cross limit, then
   stop; we don't try to refine estimate of (B0, eta2) limit. (DONE)
2. tabulated FWHMs are no longer accurate, since we have a better full model.
   Thus, recompute FWHMs when starting from best grid values (DONE)
3. when fitting w/ multiple values of epsfcn, start from best grid value
   (while keeping track of best fit found so far).  More reproducible and less
   chance that later fits get stuck in local minimum of earlier fits. (VOID)
4. Set smaller rminarc wherever possible (TBD)

(copied these over to agenda/to-dos)

See Thursday's notes about improvements to error annealing.


Thursday 2014 August 21
=======================

Summary
-------
* Only fit 1x from best grid value
* More efficient error analysis (Fitter object), testing + timing
* Skeleton code redesign

Useful: tack `.proxy-um.researchport.umd.edu` to links to get library access


Varying epsfcn systematically when fitting from best grid values
----------------------------------------------------------------
(used code commit, early morning of Friday 2014 August 22)

Previously, I let the best grid fit do the following:
1. fix best eta2, fit to find best B0 (one fit)
2. for multiple `epsfcn` values [1e-10, 1e-8, 1e-6, 1e-4], fit (eta2, B0) and
   check against previous "best" values for improvement.

Now, skip the rigmarole and run ONE fit from best grid values.  Does this still
give acceptable results?

Previously, I didn't do this because FWHMs from Fortran code didn't vary
continuously/smoothly, so nonlinear fitting went wonky and got stuck in holes.
That's what it looked like to me, but I didn't quantify/log this.

Running on SN 1006 data w/ 2 bands, results are consistent w/ Sean's Table 8
(manually checked all parameter values; chisqrs comparable or better, though
it's not an apples-to-apples comparison as I calculate FWHMs differently).

But, fitting w/ 2 data points is a bit silly since I often get chisqr == 0.
One more, let us compare single fit vs. multiple fits directly.

Filament 1, 3 bands, one fit (epsfcn = 1e-6)

    mu	    eta2	                B0	                chisqr
    0.00	24.028 (+44.637/-22.396)181.97 (+45.19/-71.29)	0.1162
    0.33	5.777 (+77.087/-4.553)	127.88 (+96.44/-26.04)	0.0568
    0.50	3.811 (+79.053/-2.713)	115.78 (+94.34/-16.76)	0.0742
    1.00	2.596 (+5.771/-1.686)	102.52 (+19.46/-9.48)	0.1415
    1.50	2.466 (+3.861/-1.367)	97.11 (+11.45/-4.54)	0.2221
    2.00	2.646 (+2.865/-1.320)	94.19 (+5.01/-2.28)	0.3021

Filament 1, 3 bands, multiple rounds of fitting (varying epsfcn)
(taken from Wednesday notes)

    mu      eta2                    B0                      chisqr
    0.00    27.998 (+72.002/-26.672)187.92 (+55.86/-79.98)  0.1149
    0.33    5.510 (+94.490/-4.490)  126.67 (+104.91/-26.29) 0.0173
    0.50    4.286 (+95.714/-3.375)  117.82 (+105.79/-20.14) 0.0399
    1.00    2.589 (+5.778/-1.679)   102.49 (+18.51/-8.67)   0.1415
    1.50    2.458 (+3.869/-1.547)   97.08 (+10.58/-5.85)    0.2221
    2.00    2.640 (+3.074/-1.620)   94.18 (+6.44/-4.34)     0.3021

Tentative result: let's ditch the multiple fits for now.  This makes the
methodology simpler / easier to explain/reproduce.  Results are good for large
mu values (B0, eta2 less correlated), and otherwise agree within error.


Test improved error annealing
-----------------------------
(used code commit, early morning of Friday 2014 August 22)

Yesterday (Weds), I refined error threshold searching (more accurate, search
farther). Today (Thurs), I refactored code using a fitter object that
stores `B0_prev` value, which speeds up fits if we are testing a monotone
sequence of eta2 values (esp. if searching outside the pre-cached grid).

### Some test results

SN 1006, Filament 1, three energy bands, mu = 0.00.
Errors for confidence interval = 68.3% (1-sigma), +1 to chi-squared value

* ONE fit from best grid values: 34 function calls (eta2=24.0283, B0=, chisqr=0.1162)
* anneal grid chisqrs, left:     43 function calls (4 eta2 values)
* brentq find threshold, left:   43 function calls
* anneal grid chisqrs, right:    32 function calls (3 eta2 values, hit grid edge)
* search past grid, right:      104 function calls (12 eta2 values)
* brentq find threshold, right:  69 function calls (>5 eta2 values)

Total function calls: 34+43+43+32+104+69 = 325.
(add ~50--100 more to marginally improve parameter fits)

SN 1006, Filament 1, 3 energy bands
left side limit: B0 = 108.24 muG; eta2 =   1.356, mu = 0.00
best fit values: B0 = 181.97 muG, eta2 =  24.028, mu = 0.00 (chisqr=0.1162)
rght side limit: B0 = 271.62 muG, eta2 = 191.251, mu = 0.00

(I have previously found better fit at B0=187.92, eta2=27.998, chisqr=0.1149)

### Observed bug (from other messing around)

Error annealing uses table FWHMs to compute chisqrs, to get limits on eta2.
But, when recalculating FWHMs with new full model code, chisqrs may change
and no longer bracket threshold as expected.  Then, spopt.brentq cannot find
root between two values and errors out.

### Conclusions?

This doubles my estimates of time needed to run fits + get errors:
~10 hours for SN 1006, 26--36 hours for Tycho's SNR.
Admittedly not too bad w/ parallel code + faster computer.

Improvements?
* grid values out to very large eta2, avoid searching outside grid
* skip `spopt.brentq`, `one_dir_root` steps and report overestimated errors.
  this requires the larger pre-tabulated grid, though.


Friday 2014 August 22
=====================

Summary
-------

Grid settings and magic numbers
-------------------------------

I changed a TON of settings to use defaults specified in SNR objects:
grid resolutions (irmax, iradmax, ixmax), cut-off energy toggle (icut), and
parameter initial guesses + hard limits for simple/full fits.
Now, `models_all.py` should be mostly free of magic constants (table generation
code aside, I suppose).

To use non-default numbers, I pass a `model_kws` argument through functions to
let fits customize grid resolution, rminarc, etc.


Error calculation methods for full model
----------------------------------------
just for fun, can we get confidence intervals out?  Vs. our grid annealing...
this assumes that lmfit approach is reasonable, which I think it's not bad but I
don't understand the details.
First run a fit to some data, then run the confidence intervals.  Compare the
results (stderr, confidence interval results):

Here: 10--23 calls to fit, 339 function calls to get confidence intervals!!!!!
(but, admittedly, it covers both B0 and eta2... so it's actually quite reasonable!)

(compare my method -- if errors are in grid, takes 80+80 function calls [l/r]
for eta2, then 80+80 for B0 = 160x2 = 320.  So time taken is about the same) 

I've got 3 different ways to estimate errors:
* lmfit confidence intervals
* `one_dir_root` without grid
* anneal fr grid, `one_dir_root`/`brentq`.  Faster but more complex

    import numpy as np
    import lmfit
    import models_all as ma
    import snr_catalog as snrcat
    snr = snrcat.make_SN1006()

    kevs = np.array([0.7, 1., 2.])
    data, eps = np.array([35.5, 31.94, 25.34]), np.array([1.73, .97, 1.71])

    # Naive fit
    res = ma.full_fit(snr, kevs, data, eps, mu=1., eta2=1., B0=150e-6,
                      eta2_free=True)
    print res.success, res.nfev, res.chisqr
    print lmfit.printfuncs.fit_report(res.params)

    # Fit based on (putative) best grid values (since grid values are outdated)
    res = ma.full_fit(snr,kevs,data,eps,mu=1.,eta2=2.857,B0=102.98e-6,
                      eta2_free=True)
    print res.success, res.nfev, res.chisqr
    print lmfit.printfuncs.fit_report(res.params)

    ci = lmfit.conf_interval(res, sigmas=(0.674,))
    print lmfit.printfuncs.report_ci(ci)
    print ci

------------------

Results (naive fit):

    True 23 0.141500732519
    [[Variables]]
         B0:       0.0001024945 +/- 4.374206e-06 (4.27%) initial =  0.00015
         eta2:     2.589774 +/- 0.9773799 (37.74%) initial =  1
         mu:       1 (fixed)
    [[Correlations]] (unreported correlations are <  0.100)
        C(B0, eta2)                  =  0.990

Results (best grid value fit)

    True 10 0.141500732512
    [[Variables]]
         B0:       0.0001024945 +/- 4.374126e-06 (4.27%) initial =  0.00010298
         eta2:     2.58977 +/- 0.9773611 (37.74%) initial =  2.857
         mu:       1 (fixed)
    [[Correlations]] (unreported correlations are <  0.100)
        C(B0, eta2)                  =  0.990

-------------------

CI from best grid value fit (but, really shouldn't matter frankly)

            67.40%     0.00%    67.40%
      B0   0.00010   0.00010   0.00011      (96.0 muG, 102.5 muG, 112.6 muG)
    eta2   1.30828   1.69929   5.28978      (? why is 0% mark different...)

    {'B0': [(0.674, 9.603169453502047e-05),
            (0.0, 0.0001024944738486107),
            (0.674, 0.00011256390784381289)],
     'eta2': [(0.674, 1.3082827379768638),
              (0.0, 1.6992912000357574),
              (0.674, 5.2897839276153755)]}


Model fitting code redesign
---------------------------

Don't spend more than 1 day's worth of work on this.
Refactoring is meant to make it easier to add next steps of functionality.

Now (re)writing error annealing for full model, in the image of my older,
more elegant method for simple model.

TODO: move presentationformatting code over somewhere...

Mostly done -- time for some nosetests, to confirm that stuff works as
expected. What are my validations?


... (added some tests... but could definitely use more)
... (enough to validate the simple fitting code works okay, though)
... (which I also verified by manual inspection)
... (now testing out new rewritten code w/ SN 1006...)



Sunday 2014 August 24
=====================

Because I have no ability to prioritize I'm still wrangling with model fit
details.  Even with the new full model code, the fits are quite inconsistent,
especially at small mu (where B0, eta2 are most strongly correlated...).

In `one_dir_root`, I change the code by not rerunning fits quite so often,
and not restoring the original (best-fit) parameter values after every
fit @ fixed (varied) eta2/B0.  The reason -- this greatly expedites the
error finding process for the full model...

I think I twiddled with some of the initial guesses, so now it's no longer
converging (test case: SN 1006, Filament 1, of course).  Ah, and of course,
I'm now searching well outside the grid.

Big problem -- as `one_dir_root` progresses -- it performs fit, then steps
forward and performs a fit.  But, often, it hasn't truly found the best fit
(and I don't fully understand why / don't know how to force it to do better).

Yeah, to be safe, I should let the fits off of grid run w/ multiple steps...

Error methods
-------------

I see a few ways to move forward with errors:

1. standard errors from numerically estimated covariance matrix
   with some knowledge of the function's behavior, and/or working from model
   equations, we could attempt to better constrain the Hessian matrix
   (i.e., use some math to get a 2nd/3rd/n-th order term expansion of what the
    true variances should be)
2. lmfit's confidence interval algorithm, ~300 function calls.
   F-test: compare null model to model w/ -1 DoF, 
    \[
        (\Delta\chi^2/\Delta DOF) / (\chi_0^2/DOF)
    \]
   * If chisqr/DOF >> 1, then we need a larger change in chi-squared to reach
     a given confidence interval... though I don't know how to show/prove this.
   * Conversely, for chisqr << 1, I think error estimates might be too small.
     Small chisqr suggests model is overfitting data, but it seems like the F
     function wouldn't capture this; rather, it would yield smaller errors...
   Looking at lmfit's code, the implementation seems quite straightforward.
   But I need to better understand this F-test (derivation, assumptions, etc.)
3. Homebrewed error searching, varying parameters to obtain an absolute
   variation in chi-squared corresponding to some confidence iterval.
   Per Numerical Recipes (Press et al.), this is only valid IF:
   * measurement errors are normally distributed
   * model is linear, or at least looks linear (within the parameter space
     enclosed by our confidence limits / uncertainties)
   NR references Lampton, Margon, & Bowyer (1976, ApJ)... it's quite relevant.
   [ADS entry](http://adsabs.harvard.edu/abs/1976ApJ...208..177L).

   Lampton, Margon, & Bowyer (1976, section IV) do test the delta chi-squared
   analysis for non-linear spectral fits using numerical simulations; for a
   given delta chi-squared contour, the confidence limits are actually kinda
   okay.  See also Margon, Lampton, Bowyer, & Cruddace (1975, ApJ) (aside: holy
   smokes those spectra... little wonder they want to be so careful about
   statistics).

   So, okay, nonlinearity is acceptable because derivations/assumptions have
   not invoked linearity (to eq'n (5) in Lampton, Margon, & Bowyer (1976))

   I suppose we can assume FWHM errors are normally distributed, for lack of
   any better expectation (that's not very compelling).
   
   But we still must deal with chisqr very large or very small.

   One ad hoc correction: the [GSL manual, section 38.11](http://www.gnu.org/
   software/gsl/manual/html_node/Example-programs-for-Nonlinear-
   Least_002dSquares-Fitting.html)
   

Mind you, no one is going to hand you a solution on a silver platter.
If there was one, it would have been written up and put in use long ago.
You must live with the fact that your errors will be at least as uncertain as
your fits (or even more so...).


(Week 13) Monday 2014 August 25
===============================

Summary
-------

Blah

Error calculation analysis
--------------------------

Studying `lmfit` source code: it uses a very similar approach to mine:
step the parameter value manually until an f-test probability is crossed,
then use `scipy.optimize.brentq` to narrow down parameter bound.

`lmfit` runs into problems because step size is too small and is invariable.
My code seems to run into problems w/ the fits converging consistently...
(I've no idea if `lmfit` would have the same problem @ that particular area of
parameter space, but I think it would... anyways, taking a page from lmfit
I have now changed the step size for `one_dir_root`-finding, testing it out now
on SN 1006 (at high mu values only)).


Oh my goodness I'm so tempted to say it's not worth it.
But I must have some numbers so that we can compare, first, and say that
scaled stderrs are tolerable.  I don't have that yet.


Parameter values bracketing error threshold are inconsistent
------------------------------------------------------------

### Example failure mode

Filament 2, mu = 1.5, 3 energy bands.  Best fit B0 = 129.6 muG.  Seek lower error limit.
* Start @ grid point B0 = 126.765 (eta2 = 0.088).
  chisqr is now above threshold, so code starts moving backwards (increase B0).
* Check grid point B0 = 127.385 (eta2 = 0.107; init = 0.088)
* Check grid point B0 = 128.886 (eta2 = 0.168; init = 0.107)
* Check grid point B0 = 129.739 (eta2 = 0.333; init = 0.168)
  We found B0 below chisqr threshold (but passed best fit B0, FIX this!)
* Call brentq rootfinder:
  * Check limit pt B0 = 128.886 (eta2 = 0.209; init = 0.333) (changed!!!)
  * Check limit pt B0 = 129.739 (eta2 = 0.333; init = 0.209) (same)

`sp.optimize.brentq` throws a ValueError, as the limits now have same sign.
Here, even if I searched between the found edge (128.886) and the best fit B0
(129.6), I would have had the same result.

Made some ad hoc changes, drawing from `lmfit` approach.
* If error search doesn't move on grid, just search btwn grid crossing and best fit
* In both `one_dir_root` and `anneal`, if brentq gives ValueError
Earlier I changed step size eps to follow lmfit.
Now I added try/except block to catch when the fits differ, and make a feeble
effort to correct for it...

A couple more heuristic tweaks to step sizing etc... to be suitable for my
usages... really ad hoc :/


Mysterious resolution errors (sizable B0/eta2, mu=0)
----------------------------------------------------
See debugging log for information on this. Bug identified 2014 August 25,
but not addressed until a bit later.


Bad initial guesses kill rootfinding
------------------------------------

Searching outside pre-tabulated grid requires new fits, initialized from
previous best fit parameters.  But, once we find parameter values bracketing
our error threshold, we must often make a large jump in parameter space.

Example (filament 1, mu=0.33): best fit is B0=127.45 muG, eta2=5.678
Bounding B0 from below was easy (B0 = 100.454 muG lower limit),
but bounding B0 from above fails.
We hit grid edge with B0 = 231.068 muG, eta2=99.018; we search until
point B0 = 523.554, eta2=16641.  Then, we search between grid edge and new
bound (from `one_dir_root`).

But, starting fit with B0=231.068 muG, eta2=16641 fails horribly!

The obvious contrast is to search between two adjacent search points, 1 which
was found to be below threshold, 1 which was found above.
But that gave a similar problem, where the fits would change and the points
would no longer bracket correctly.

It's kind of depressing to stare at your rootfinder marching towards an
untimely demise.

Ad hoc solutions:
* go back to searching btwn two nearest points possible.
  when we get a ValueError, push the two points out farther left/right.
* Keep the best fit parameters at each point in time.  Reinitialize
  the fitter w/ these values whenever we must jump backwards.
  (when feeding numbers to brentq, feed numbers for 1st limit pt
  then PRAY that brentq can find the right fit on 2nd limit pt...)


Modifying `lmfit.conf_interval` for our use
-------------------------------------------
Possible.

1. `prob_func` is easy to implement (I don't know if d(chisqr) approach should
   be favored or disfavored vice f-test)
2. really need adaptive step size for my model function
   (lmfit confidence interval gums up too often at 200 iterations)
3. in `lmfit.confidence.find_limit(...)`, we have to search farther away from
   limit if brentq raises a ValueError (lines 230--245)

Otherwise, this does everything I've been laboriously reimplementing.

If I have free time (ha. ha. ha.), I could update `lmfit`'s `confidence.py`
and submit a pull request.  Would require:
* add smart stepping functionality
* run `lmfit`'s test suite
* review the f-test statistics
Anyways I forked it for now.

Also check up on that default sigma value, for f-test.  Does it correspond
to 1-sigma?  (0.674 vs. 0.683)  Given defaults for 2,3-sigma are 0.95, 0.997;
these DO correspond to the 2,3-sigma cases of the normal distribution.


Relevant literature on statistical tests
----------------------------------------
Because I'm accumulating a pile of browser links which I haven't fully parsed.
I probably will never get to parse everything.  But at least I can be more aware
of what I don't know.

* W. H. Press himself on nonlinear fitting:
  [link](http://nr.com/CS395T/lectures2009/2009_10_NonlinearLeastSquares.pdf)
* "On the unreliability of fitting":
  [link](http://hea-www.cfa.harvard.edu/astrostat/slog/groundtruth.info/
  AstroStat/slog/2007/on-the-unreliability-of-fitting/index.html)
* Related post on bias in fitting:
  [link](http://hea-www.cfa.harvard.edu/astrostat/slog/groundtruth.info/
  AstroStat/slog/2007/astroph-07054199/index.html)
* Lampton, Margon, Bowyer (1976), ApJ
  [ADS](http://adsabs.harvard.edu/abs/1976ApJ...208..177L)
* Protassov, van Dyk, Connors, Kashyap, & Siemiginowska (2002), ApJ
  [ADS](http://adsabs.harvard.edu/abs/2002ApJ...571..545P).
  I like that the running head is "Fallible f-test".
* stats.SE on degrees of freedom:
  [link](http://stats.stackexchange.com/a/17148)



Tuesday, 2014 August 26
=======================

Summary
-------
* Review host of issues, bugs, notes from yesterday.
* Check lmfit standard error scaling (cov. matrix scaling in scipy/etc)
* Refactor / clean up formatting code in table/plot ipynb files


Verify / check covariance matrix scaling
----------------------------------------

If we're going to consider using fit standard errors (i.e., sqrt of diag of
covariance matrix), we need to know how it's calculated.
Scipy tells us:
* "cov\_x is a Jacobian approximation to the Hessian of the least squares
   objective function."
* "This matrix must be multiplied by the residual variance to get the covariance
   of the parameter estimates – see curve\_fit"

The distilled, simplified code follows.
In short, `curve_fit` and `lmfit` both scale the covariance matrix by
reduced chisqr already.  In fact, I may need to remove this scaling... ?!
Where the GNU software library recommends scaling by a factor of redchisqr,
the example does appear to be a weighted fit already, without preemptive scaling
(of `lmfit`/`curve_fit`) applied.

Possible courses of action:
* remove scaling from lmfit's reported values.
  parameter errors ASSUME data errors = 1-sigma, but consider data errors
* leave scaling in place.
  parameter errors are, scaled as if data errors gave chisqr=1
* remove scaling, and replace with max(1, sqrt(redchi)) (see GNU sci. lib.)
  I don't konw how the hell to justify this, really.

ALSO (aside): for cov matrix to be valid, errors must be 1-sigma errors on input
data points.  Dammmmnit.  I'll leave it and make a note for now.


`leastsq` (SciPy 0.14.0):

    perm = take(eye(n), ipvt - 1, 0)
    r = triu(transpose(fjac)[:n, :])
    R = dot(r, perm)
    cov_x = inv(dot(transpose(R), R))

`curve_fit` (SciPy 0.14.0):
    
    # cov_x = scipy.optimize.leastsq output
    # func = leastsq objective function, returns residuals
    if not absolute_sigma:
        s_sq = (asarray(func(popt,*args))**2).sum() / (len(ydata)-len(p0))
        pcov = pcov * s_sq

`leastsq` (lmfit 0.8.0-rc2):
    
    fjac = transpose(fjac) / take(grad, ipvt - 1)  # What is this for?
    
    perm = take(eye(n), ipvt - 1, 0)
    r = triu(fjac[:n, :])
    rvec = dot(r, perm)
    covar = inv(dot(transpose(rvec), rvec))
    
    sum_sqr = (resid**2).sum()
    nfree = ndata - nvarys
    covar = covar * sum_sqr / nfree

GNU software library:

    covar = inv(J^T J) by QR decomposition (if A = QR, A^{-1} = R^{-1}Q.T)
    
    If minimisation is weighted, covar gives correct stat. weights
    If minimisation is unweighted, covar should be multiplied by ssqrred/dof



Wednesday 2014 August 27
========================
(needs cleanup)

Summary
-------


Code review
-----------
How many layers of abstraction have I built up?
1. full model code, in Python and Fortran
2. `models.py` and `snr_catalog.py` wrap the fitting model / SNR details
3. `models_exec.py` calls the fits, builds SNRs
4. `models_disp.py` factory function for fit generation, and tables/plots
5. `*.ipynb` actually call functions with various data, settings, etc

Welp.  Well, it works nicely now, and I think the documentation is pretty clear
and robust.  So it's easier to jump up/down levels of abstraction.

Here I `git commit` before moving to another debugging round (trying hard
to make error calculation more robust).


Debug "extreme-eta2-glitchiness"
--------------------------------
See Monday August 25 for notes on identification


More robust error finding
-------------------------
This should partially sidestep the extreme values bug...
but, no matter what, it will be a killer if encountered.

Refactored code; brentq now searches only between two closest points at first.
Then, if fails, it tries again with backup points.
Then, if fails, it prints a warning, returns ZERO error, and moves on.

1. if fitting between two GRID points

....
long pause of coding 

Did a lot of grid point twiddling.
The approach is now to use the smallest bracketing interval possible


So, I changed the bracketing behavior of brentq...
Simple fit is using initial step size None (defaults to 1% of dist between
initial and limit x).
Full fit, following lmfit somewhat, uses the parameter stderr to start.


Thursday 2014 August 28
=======================

Summary
-------
* Picked up GSFC badge
* Post mortem of overnight attempt at error computation
* Parallelizing iPython (led to even cleaner code, now!)

Note: when possible, try to clear easily-regenerated iPython output cells
I will still keep longer outputs (>10 minutes or so), which may be worth having
backups of anyways.


Overnight error annealing run
-----------------------------
Last night I tried running SN 1006 computations.  An incomparable trainwreck
ensued.  I went through verbose output and took notes for next round of
debugging (if necessary...) (notes not included in today's commit).

iPython parallelization
-----------------------
With some twiddling, code runs nicely in parallel (spent a few hours hacking
around with this).  Need `pip install dill`, mind you.

DirectView's `map` (or, `map_sync`/`map_async`) does not play nice with inline
matplotlib plots.  Compare to line magic `%px` which manages to make it work.
I don't really understand why, and it's a bit much to dig into.
The best I can do is to display plots in separate windows.

Workaround -- don't plot in parallel code.  Save the numbers as output, and
plot after the fact


Friday 2014 August 29
=====================

Summary
-------
* ipynb to check resolution errors
* Re-lowered eta2 upper limit (now: eta2 max = 1e5, B0 max = 1e-2)

Remark: the 'I' in IPython is supposed to be capitalized.  Woops.

Meeting
-------
* Write (you should have started weeks ago, or done it meanwhile)
  Send material anytime (goal: outline + methods within a week)
* Perhaps, ask Keith Arnaud about errors -- low priority...
* Push off magnetic damping until later (address when we touch on it in
  discussion)

More friggin debugging
----------------------
Finished code post-mortem.  Fixed some egregious bugs in error root finding.

Resolution analysis
-------------------
Can't put it off any longer, because I need new tables for error computation.

Happily, with some refactoring while coding (and some thought to how it will be
used/reused)... it wasn't so bad.  Helpful to think of the iPython notebook as
just that, a notebook.

Observation: the largest errors are associated with highest energy bands
and negligible diffusion (eta2 = 1e-16), interesting.  Why?
Hypothesis: we're probably not resolving FWHMs as well (high energy, low
diffusion --> narrowest possible FWHMs).

In Tycho, this is a problem.  Default rminarc of 20, when some of the smallest
FWHMs are about 1-2 arcsec, drives the errors way up.

Possible solutions/improvements, for Tycho's SNR:
* double Tycho iradmax to 200.  Set rminarc default to tighter values,
  e.g., `[20., 16.5, 13.5, 14., 14.]`.
* adaptively tighten `rminarc` in `models.py`.  Leave `iradmax = 100` and
  `rminarc=20` (now rminarc doesn't matter so much)

Adaptive rminarc
----------------

Updated `models.py`.  Procedure:
* compute FWHMs with specifid `rminarc`, but use 1/2 `iradmax` resolution
* set `rminarc = fwhms * fudge_factor` (I take 1.2 as default)
* recompute FWHMs with new `rminarc` at full `iradmax` resolution

The idea is that, compared to e.g., doubling `iradmax` naively, this is
slightly faster and achieves better precision.
Side benefit that it may be more robust too!

The code will be slower than before, but I hope it's worthwhile.
We couldn't do this before with the strange Fortran floating point (?) bug.
But, now we can!

Testing on the resolution analysis code, drops the MAX resolution error from
initially ~10-30% (iradmax=100, rminarc=20) to 0.1-0.5%, (iradmax=100).
This is ~2 orders of magnitude improvement, but not even doubling the
computation time (a pinch below).  Hallelujah!


New tables
----------




