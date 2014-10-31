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

* Week 10 - (8/4) full model grid-best-"fits" (test fits, error annealing w/
            SN 1006). Port full model code to Python, optimize for speed
* Week 11 - (8/11) one week break

* Week 12 - (8/18) LaTeX table code; more precise manual error calculations.
            Debug/test new full model
* Week 13 - (8/25) refactor model exec/disp code.  Debug error calculations
            extensively.  Back on GSFC campus.
* Week 14 - (9/1) run suite of error calculations; start paper text/outlining

* Week 15 - (9/8) manuscript tables, varied FWHM calculations,
            first Kepler picks; paper text
* Week 16 - (9/15) flesh out paper, many new plots/tables (variant fwhms,etc).
            New Tycho region selections
* Week 17 - (9/22) paper writing, cleanup/data for Tycho regions-5
            Cull Kepler regions, trial Cas A regions
* Week 18 - (9/29) AAS abstract. Trim paper, flesh out results/disc/conclusion.
            Port B-damping code.
* Week 19 - (10/6) Investigate srcut break-D relation; prelim B-damping fits.
            Check code resolutions, apply new settings.  Make B-damping table
* Week 20 - (10/13) B-damping fits (Tycho, SN 1006).  e- distr derivations?



* Week 21 - (10/20) Emails w/ Sean, paper writeup on damping results
* Week 22 - (10/27) Tycho regions-6, fits w/ srcutlog eta2



* Week 23 - (11/3)
* Week 24 - (11/10)
* Week 25 - (11/17)
* Week 26 - (11/24)
* Week 27 - (12/1)
* Week 28 - (12/8)
* Week 29 - (12/15)
* Week 30 - (12/22, 12/23 only) 

* Week 31 - (12/29)
* Week 32 - (1/5) AAS winter meeting

(week 10 included for continuity)


(Week 11) Monday 2014 August 11 -- Friday 2014 August 15
========================================================

On vacation


(Week 12) Monday 2014 August 18
===============================

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

    mu      eta2                    B0                  chisqr
    0.00    24.028 (+44.637/-22.396)181.97 (+45.19/-71.29)  0.1162
    0.33    5.777 (+77.087/-4.553)  127.88 (+96.44/-26.04)  0.0568
    0.50    3.811 (+79.053/-2.713)  115.78 (+94.34/-16.76)  0.0742
    1.00    2.596 (+5.771/-1.686)   102.52 (+19.46/-9.48)   0.1415
    1.50    2.466 (+3.861/-1.367)   97.11 (+11.45/-4.54)    0.2221
    2.00    2.646 (+2.865/-1.320)   94.19 (+5.01/-2.28) 0.3021

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
* Fixed misc. error finding bugs; FWHM calculation uses bisect over brentq now
* Adaptive rminarc calculation in full model code

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

Extreme value bug
-----------------
See log.  Solution is to replace `scipy.optimize.brentq` with
`scipy.optimize.bisect` for rootfinding.  However, the code easily goes bonkers
at such extreme values (intensity values rapidly hit the limit on representable
numbers...).  So, better to limit parameter values appropriately, as much as
possible...

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

We could move to ipynb to take advantage of notebook interface, esp. for
parallelization (benefit from, e.g., load balancing would be fairly marginal).
Main benefit would be the nicer interactive interface.

But, current code/procedure works all right.  It ain't broke, don't fix it.
I cleaned/updated the code slightly, but it's mostly unchanged.
It seems to work like a charm.

### Shock velocites to tabulate

Previously used shock velocities:
4.21e8, 4.52e8, 4.83e8, 5.14e8

Min, max shock velocities from regions-4 assuming d = 3 kpc:
4.50e8, 5.20e8

Today, let's use: [4.59e8, 4.76e8, 4.94e8, 5.11e8] which should give roughly
even spacing for points in [4.50e8, 5.20e8].  Generate with code snippet:

    f = np.array([1., 3., 5., 7.]) / 8
    v = 4500 + f * (5200 - 4500)

This should give slightly closer spacing overall; of previously generate
tables, the lowest shock veloc. was never used by my regions.

### Gridding parameters, time

Currently covering 4 decades in eta2 ( the log spacing) w/ 50 points.
Try covering 5 decades, -2 to 3, w/ let's say 60 points.
Testing out eta2 = 1e3 on Tycho, looks okay, I'll go ahead and do it.

* 6 mu, 110 eta2, 20+ B0
* 4 shock velocity values
* 4 sec/call  (being pretty conservative)
 = `6*110*40*4*4 / 3` (parallel) = 140800 sec = 39.1 hrs (!)

Expect to be done Sunday evening (8:30pm ish), possibly earlier since I think
this should be a conservative time estimate.  All right, here goes.
Once tables are generated, merge them together and let error annealing run
overnight.  Do the same for SN 1006.

### Table generation procedure (2014 August 30 tables)

Grid generation procedure has changed since 2014 July 30.
New procedure looks like (adapted from `rsch-notes.md`):
1. Write script to call table generation method w/ desired parameters
2. Check activity monitor for processes hogging CPU, etc.
   Change energy settings so that computer doesn't sleep (display off only).
3. Call code as:

    python fullmodel_recompile.py
    git add -A
    git commit
    python one_time_script.py  # whatever you want to call it

Note: initial guess for eta2 always fails badly on eta2 = 0, or eta2 extremely
small (`width_dump` function goes singular?).  Not to worry, code will adapt.

4. Search log files for errors (ignore caps)
5. Move tables and log files to /tables and set as read-only
6. Merge tables if needed (can do this manually, or whatever)

Start: 2014 August 30, 05:30am
Finish: ???


(Week 14) Monday 2014 September 1
=================================

(Labor Day)

Tables finished ~7pm, just about as expected.

Tested out Tycho fits from new grid, looked okay. Started error annealing
calculation running overnight.

For some reason, after ~7 hrs one of the cores stopped.  Weird behavior going
on -- mu values are jumping around in the output log?!  I'm not sure what's up.
Wait until all cores done and look at display output, tables to see...


Tuesday 2014 September 2
========================

Summary
-------

* New SN 1006 grid, Tycho error calculations running (start-to-finish trial)
* Start some paper text

Misc: installed Pathogen, ConqueTerm for MacVim usage.  iTerm2 or something
else probably makes more sense, but this works for now.
Cleaned misc. files, repository structure. General housekeeping.


New tables, calculations
------------------------

Start generating new SN 1006 grid around 11am. `rminarc`-too-small errors are
popping up at the preliminary step (calculation w/ small `iradmax`), but the
output FWHMs look okay (no weird numbers).

Expected time? 6 x 60 x 40 x 3 sec/call = 12 hours (not running in parallel).
Still running Tycho error calculations so making use of 3 available CPUs.

Tycho errors 

Quick tally: region 1, mu=0.
* Fit from grid:        20 function calls
* B0, lower bound:      186 function calls (moved 1x on grid)
* B0, upper bound:      125 function calls (moved 1x backwards on grid)
* eta2, lower bnd:      150 function calls (moved 1x on grid)
* eta2, upper bnd:      181 function calls (moved 1x backwards on grid)
* Total:                662 function calls = 55.17 minutes (@ 5 seconds/call)

* Best fit: 784.180 muG, eta2=13.664
* B0 lower bound: 711.870 muG (eta2=12.051)
* B0 upper bound: 908.343 muG (eta2=38.075)
* eta2 lower bnd: 12.023 (B0 = 720.636)
* eta2 upper bnd: 38.253 (B0 = 916.589)

Time: 13 regions x 6 mu values x 662 calls x 5 sec/call = 258180 sec,
divide by 3 for parallel set-up (maybe, 2 since some cores stopped/died?!)
Expect: 24 to 36 hrs (!).  Holy smokes.  So it should be done by tomorrow
morning...


Manuscript
----------

LaTeX/BibTeX setup ([relevant](http://www.astrobetter.com/install-latex)).
BibDesk, latexdiff are included w/MacTeX.  
* AASTeX, emulateapj, `apj.bst` (fr. AstroNat).
* J. Sick's adsbibdesk (3.2.0); Automator for OS X + command line
* BibDesk, latexdiff included w/ MacTeX.

I use cite key format `%a1%Y%u0` (e.g., `\citep{ressler2014}`).
Disable pdf downloads with adsbibdesk (in `~/.adsbibdesk`).
Do not link files to BibTeX file -- messes up plaintext output and further
clutters git diffs.
Imported most (i.e., relevant) references from Sean's paper

Auto-generated BibTeX files are going to be ugly as sin -- obviously, will be
uglier than doing by hand.  Accept the trade-off: no need to prettify when it
is parsed fine by BibTeX and looks okay in final paper.

Utilities for collaborative editing:
* trackchanges (F. Salfner, [sourceforge](http://trackchanges.sourceforge.net))
* commands for collaborators to markup document
  ([ex.](http://emckiernan.wordpress.com/2012/11/28/tracking-changes-in-latex))
* latexdiff for indiv. revisions (I keep track)
* ShareLaTeX (10 collaborators @ $6/mo, would need to get people on board)

Tentative workflow -- work directly in `paper.tex`, no jumping between files.
Use adsbibdesk, BibDesk to interface with `.bib` (avoid temptation to manually
twiddle).  Using `trackchanges.sty` by F. Salfner.

Best practices:
* commit @ every iteration emailed out for comments (mention in commit)
* commit after all comments are added


Wednesday 2014 September 03
===========================

Summary
-------
* Tycho error annealing finished, investigated bugs
* Start SN 1006 error annealing
* Paper outlining


Tycho error annealing post-mortem
---------------------------------

Finished, though I neglected to record the time it finished.  Doesn't matter
since we hit an error in the last engine (index 3, of 0--3).

Added annotation to help in debugging... (determine why the index error is
occurring).  Currently, recomputing errors for region 11 to try to trace error.

I'm perplexed because I genuinely don't expect an IndexError where the
exception popped up.

From another calculation (run through today) -- it looks like the error finding
just dies?!  But everything else behaves...

SN 1006 new grid, full model/error calculation
----------------------------------------------

I count 1019 logged ERRORs in code, all identical -- correspond to the initial
estimated calculation that allows us to refine our guess for `rminarc`.
No instances of "Error!", which would indicate that the 2nd FWHM computation
(the final result) went awry.

So, new SN1006 table looks good (`sn1006_gen_2014-09-02_grid_6-110-20.pkl`).
Start error annealing, allow to run overnight (3 bands first).

Paper outlining
---------------
Finished infrastructure / set-up / messing around with AASTeX.
Started sketching out introduction, methods, moving text over from
poster/notes.


Thursday 2014 September 4
=========================

Summary
-------
* More paper writing/outlining (messing w/ tables, text, etc)
* SN 1006 tables (started yesterday) failed on error (edge case)
  I believe `check_f_init` wasn't initialized.

Friday 2014 September 5
=======================

Summary
-------
* Disc. w/ Brian -- move forward with numbers, let errors come later
  For Kepler, stick to the south.  Writing -- just mark up with red pen.
* Spent a little time on advective solution to transport eq'n.  Maybe
  important, esp. as large (resolution) errors seem to appear in the advective
  case.


Saturday 2014 September 6
=========================

Following Rob's email -- quick calculation on Tycho to see if rim positions
have any energy dependence in our data.  Appears not.


(Week 15) Monday 2014 September 8
=================================

JVGR...

Tuesday 2014 September 9
========================

Summary
-------
* Rewrite error annealing code, much simpler now...
* Table formatting galore

More error annealing debugging
------------------------------

Annealing -- no longer use parameter stderr as initial step size.  Too large,
often.  Let the rootfinder start w/ small guess, and adapt step.
(observed in Tycho, at edge cases, steps too large)

Annealing -- stop with the grid stepping.  Now, only use the grid with chi-sqr
values to get the place to try and start searching.

Completely revamped and simplified annealing now.  No longer actually anneals.

SN 1006 output
--------------

Went through and pulled out relevant numbers.  Data storage organization
pending...


SN 1006, two band fitting
-------------------------

Now generating new tables/errors for SN 1006 with two bands, using new error
computation code.
How long did this take?

    [u'ok', u'ok', u'ok', u'ok']
    16347.85151
    4:26:20.616579
    2:43:13.325094
    2:31:40.705533
    2:17:11.627385

How many function calls were needed?  Compare to my old approach.
Filament 1, mu = 0, 2 energy bands

* fit from grid: 48 calls
* lower bound on B0: 205 calls (18+15 to move 1x from grid, 172 bracket)
* upper bound on B0: 132 calls (18+18 to move 1x from grid, 96 to bracket)
* lower bound, eta2: 163 calls (16+13, move 1x from grid; 134 to bracket)
* upper bound, eta2: 245 calls (14+11, move 1x from grid; 220 to bracket)

Total calls: 48 + 205 + 132 + 163 + 245 = 793 calls

At ~2 sec/call, 6 mu values; 1 filament takes ~9500 sec = 2.64 hr.
In fact it's a bit faster -- maybe fits at higher mu are easier.
But the estimate is in the right ballpark.
Compared to my previous estimate, of 662 calls/filament/mu (Tycho), using 662
calls here actually gives about the right amount of time (at 2 seconds/call).

Conclusion: new error annealing (rather, searching instead of annealing) is
likely the same speed as before, and seems about as good -- and is much less
error prone.


Tycho, fitting filament-averaged FWHMs
--------------------------------------

Time (stderr only):

    [u'ok', u'ok', u'ok', u'ok']
    1660.824469
    0:23:50.932605
    0:12:00.476653
    0:13:22.696450
    0:07:46.230326

Skimmed stdout log, looks good (I think) but I didn't look very carefully.
Mainly no explosions of errors... just a few while fitting.


Typesetting nice tables
-----------------------

Dear goodness this is a bit tedious

Vim snippets for manual table spacing (this is _bad_ practice).
Adds more space between best fit value and super/sub-scripted errors:

    %s/\^{+/^{\\,+/g
    %s/_{-/_{\\,-/g

Okay.  See committed version of paper for a bunch of tables.
Emailed first round draft to Brian/Rob


Wednesday 2014 September 10
===========================

Summary
-------
* ...

Tycho, average filaments
------------------------
Started error computation around 6:05 am
expect it will take ~4 hours... sigh.  x2 because 4 engines for 5 filaments, so
one engine is gonna work overtime.  

(edit: 10.5 hours?! dammit.  Well, close enough to my estimate of 8 hours.  But
probably closer to 5-6 hrs/filament it seems)

Misc things
-----------
* Started playing with Kepler, picking very fine regions -- looks like we can get
  some good data, I think.
* Considering error calculation, presentation for FWHMs...
  my plan is to use a geometric average, which preserves energy dependence.
  estimates of B field, etc may be less meaningful?  Not certain.
* Error on geometric average is tricky.
* Short UMD trip


Thursday 2014 September 11
==========================

Summary
-------
* Kepler and JVGR...
* Misc paper writing

Kepler region picking
---------------------
This is relaxing work...

Quick thoughts:
* could we smooth images to get more robust rims -- reduce count fluctuation?
* select one region, then generate a whole bunch with artificially varied pixel
  offsets, rotation angles -- then we get an ensemble of profile fits and FWHMs

Finished selecting regions


Friday 2014 September 12
========================

Summary
-------
* Consolidate / add more table outputs
* Averaging FWHMs and getting statistically meaningful errors
* Cleaned agenda to-dos (why is there so much backlogged crud on here)
* Start running geometric FWHM + error calculations

Tycho fits with eta2=1 fixed
----------------------------

Updated / patched code to work with either eta2 or B0 fixed.
Now running on Tycho averaged filaments, for paper draft presentation

Realization -- does it make sense to give errors when there is only one free
parameter?  Not really... so don't waste time computing error, actually.
Because the results won't really make sense (and they let eta2 run free cuz I
haven't fixed the code).  Made note in code docstring

Some engines appear stuck/hung up on error finding.  Terminating and re-running
without error calculations...

Exercise (sanity check of error fit)
------------------------------------

Try plugging numbers into full model, manually, on Tycho Filament 1 data
(arithmetic average of FWHMs) with mu=2, using recent output data.

Assuming the given best fit is indeed the best, the numbers look pretty good.


Averaging FWHMs and computing errors
------------------------------------

I discussed this over with Brian a few times but I'm really uncertain as how to
move forward -- this is where we should throw it in the paper, so that we can
discuss with Steve/Sean better...

Typed up short justification for geometric error
Added explanation of error calculation to tex file.
(edit: removed Sept 16)

Determined how to get error on geometric means (after some fumbling, staring at
older papers, etc), as follows:
for small sample sizes, mean of n samples follows (n-1) student's t-distr.
Using this distr, we get new 1,2,3-sigma errors from stderr by adding a
multiplicative factor. (these correspond to locations of %-spread on
t-distribution, but are farther out because of the small sample size).
In fact, this is important/relevant for our arithmetic means too.

Geometric FWHM + error calculations
-----------------------------------
Git commit settings, started running tonight at 01:00a, 2014 Saturday Sept. 13.


Saturday-Sunday 2014 September 13/14
====================================

Summary
-------
* Work on Kepler pipeline

Kepler data work
----------------

Generated counts files for Kepler (from `merged_evt.fits`).
Working off of `regions-1/regions-1.reg`.

Now, just retracing pipeline.  Thank goodness for documentation.
* Generate profile data, for mosaics and counts files:
* Generate background spectra.  Currently twiddling stuff by hand:
  set the ObsIDs, set the weights, etc...  Might have to wrestle with it a bit.
* Copy and modify profile fitting notebook Kepler.
  Made several ad hoc tweaks to code, Kepler notebook to get it to work
  (Tycho notebook no longer functional, but easy to fix).
* Generate region spectra (w/ upstream/downstream cuts)
  running now

WARNING: CURRENT FWHM PKL FOR KEPLER REGIONS-1 DATA IS TWIDDLED
ONLY TEMPORARY -- DO NOT USE FOR DERIVED DATA, MODELING, ETC

I changed data in `regdict[n][lab]['meas']['fwhm']` (and limits/errors) to use
the `fwhmc` data -- i.e., same profile fit, but taking the maximum at the
highest data point instead of highest fit point.
Reason -- with default fit settings, one of the regions couldn't find a FWHM in
any bands (region 13, best fit strongly /undershoots/ peak).

### Sidebar on coordinate inconsistency

I think I've been using inconsistent coordinate systems etc... background
links, spectra, regions, etc may all be slightly off (physical coords may
differ w/ single obsID vs. merged obsID...).  Issue identified w/ pipeline from
back in July or so.

Brian's `merge_spec` scripts do:
* convert wcs to physical CIAO coords for each ObsID using the tycho imgs.
* run specextract from phys CIAO coords, on EACH ObsID's `_repro_evt2.fits`
* merge spectra together w/ `ftools` stuff

Each ObsID has a different physical coords -- WCS transformation.

My region files are WCS (fk5).  When I convert to Chandra physical coordinates,
I'm using the `merged_evt.fits` physical coordinates. But, feeding these
physical coordinates into specextract on a single ObsID will give spectra from
the wrong location!  I believe this was the issue w/ my background spectra,
taken from only one ObsID; the regular spectra were merged from the 750ks
exposure.

Maybe there was another issue, who knows... but I think that was the big one.

Anyways, point is -- main solution is to just run the scripts to get spectra
from ALL ObsIDs, for consistency.  We'd have to do it for publication figures
anyways, and we're not running this that frequently.

Model fit outputs
-----------------

Ah, deferring for now... need ad hoc tables first, to have numbers assembled.
Added very ad hoc code to generate FWHM tables.

But, I must revisit this...
Not committing for now.  Will revise code tomorrow.


(Week 16) Monday 2014 September 15
==================================

Summary
-------
* Morning meeting - how to best present results, organize text, interpret
  numbers...

Meeting
-------

* Lots of comments thanks to Brian
* Spectra stuff before FWHMs.  Present tables with chi-squares from fits
  before/after filament.  Did we come to a consensus on how to do the
  downstream fits?  Maybe not.
* Simple model -- present results for 1 filament for comparison, then
  give sentences -- it doesn't make a difference, so might as well use the full 
  result.
* Reporting model fit chi-squares and plots -- the chi-square is so large as to
  be practically meaningless (it doesn't tell us whether our model even makes
  sense or not, we're just trying to make the shoe fit).  Could we focus on
  just the plots of best fits, and derived errors?
* What to do about tables, errors -- move forward by averaging within filament.
* Magnetic field -- we have strong fields, no matter what we do.

* Kepler -- focus on one good filament.  Show some other regions for
  comparison, see how they are similar/differ.

But, paper is absolutely first priority.  Get a more complete version to
Rob/Brian asap (Rob refrained from commenting since it was so un-fleshed out at
this stage), then we can send around to collaborators.

Paper writing
-------------

Added more flesh to introduction, spectra, filled in some number of citations
and details.

FWHM data organization
----------------------

Come back here.  I think I have to be more data-centric, less pipeline
focused.  Data, numbers are more permanent than may be thought... if only by
pure inertia (that once generated, there is a delay to generate new data).

Plan: figure out how to save the data smartly.  Manually convert some old data
over (and make a note of it, and how to regenerate it).

Then, write/modify scripts to generate publication-quality figures.
These should be separated from the data processing pipeline -- but can be
invoked by the pipeline.


Cas A set-up
------------
Start running `reproject_obs` as I write.


Tuesday 2014 September 16
=========================

Summary
-------
* Paper writing - fleshing out + layout planning

Paper writing
-------------
Finish incorporating Brian's comments, ideas from mtg yesterday into text.
Fleshed out more of introduction, methods, results.  Cleaned, restructured
layout in anticipation of planned tables/numbers to report.

Plan is to go back to code, now, with better idea of what must be generated.
Then, create all necessary tables and SAVE the data.

Wednesday 2014 September 17
===========================

Summary
-------
* Generate more Cas A, Kepler energy band images
* FWHM code rewriting

Cas A energy band images
------------------------

Planning to run CIAO `flux_obs` on reprojected Cas A data, but something is
amiss -- the edges of the remnant are cut off in the reprojected files,
and in the repro'd evt2 files.  But, after testing `dmcopy` with a specified xy
grid, it looks okay (I hope).

Runnning CIAO `flux_obs` on reprojected Cas A data, then `dmcopy`.

Code organization
-----------------

Pulled out script to process profiles and identify fit domains.
Massively restructured FWHM processing and fitting -- mainly threw out
massively cluttered region dictionary in favor of smaller files with less
clutter.

Pull FWHM table generating code to new folder for publication tables/figures.
Generate new, full table w/ new region dictionaries; added to paper.
Includes global averages of FWHMs and m\_E values, though std. deviations are
quite large...

Added config log method to SupernovaRemnant objects, to spew out all internal
parameters indiscriminately.

Saving best fit parameters and errors from full model fits.


Thursday 2014 September 18
==========================

Summary
-------
* Code to save/read simple/full model fitting outputs to files, including
  metadata
* Regenerate profile/spectra figures for paper; generate spectrum fit table

Code organization II
--------------------

Assembled data serialization for our model fit outputs.
Code to generate tables, merge tables for model fit parameters

Rewrote LatexTable class -- now much simpler, and can merge tables together (so
we can stack regions side by side).

Finally, what's the best way to handle plot/table generation?
Maybe generate new scripts for each dataset, and copy-paste as needed.


Figures/tables
--------------

Reviewed and resized profile plots, spectra plots for paper.
Added (almost) paper-ready plots, spectrum fit table.

Reviewing region selections: Regions 9, 11 should not be blacklisted.
Using FWHM code that 1. does not subtract background, 2. does not cap FWHM
calculation (so, most likely to accept these regions and find a FWHM).

Although the peaks are iffy and the FWHMs are close to the limits of detection,
we can't throw these data points away.  Paper text and code blacklist updated.
Note also that the manually specified blacklist may not be all inclusive.

Minor tweaks/text changes to paper.
Assembled average of full model best fits by hand.


Friday 2014 September 19
========================

Summary
-------
* Tables -- assemble some manually, added code for model fit tables.
  I think, good enough to send around

Table scripting
---------------
Assembled tables of individual region fits by hand, using new output pkls.
Added code to generate tables, + changed last column to reduced chi-squared.
Added necessary code to speed up table generation and merging.

Added SNR to pickled model fit output, to enable subsequent work to recreate
the fit exactly (mainly for plots).

Talked w/ Brian about a few paper details, and on balancing coding vs manual
work.  As discussed yesterday -- we don't need to create a smooth continuous
pipeline here.  Better to just do the work in chunk; move things manually if
needed.  We only need to run this a few times.  Added some small corrections
noted by Brian.


Tycho FWHMS with manual cap + background subtraction
----------------------------------------------------

Now running full model fits on these guys, for comparison.
Goal is just to get chisqrs out.  Quickly generate a table for all regions and
compare to our current table.


Sunday 2014 September 21
========================

Summary
-------
* Add tables w/ modified FWHM calculation, simple model fits
* Paper text (intro, transport models, results)
* Select regions-5 (see notes), run through pipeline

The new regions look workable.  Bad, perhaps, but the level of badness is
comparable to that of the original regions anyways...


(Week 17) Monday 2014 Sepetember 22
===================================

Summary
-------
* Generate tables for new region-5 FWHMs, additional tables/numbers for
  modified FWHM (cap + background subtraction) calculation, regions 4 and 5.
* Review paper layout for intro / transport models
* Meeting -- main point of consternation (call in the theorists) is, well, eta2

After generating slew of tables, started walking through intro/transport text.


Meeting remarks
---------------
For myself: aim to send material that can be distributed to Sean/Steve, by next
Monday... (earlier if possible... so Brian/Rob can look it over by Monday).

Rob's idea: to get an idea of the sensitivity of our fitting procedure, our
modeling procedure in general, run fits on fake data.
(adding noise is a separate thing).

Ideas of slicing up regions (lots of overlapping things) -- Rob suggested, move
along and find most homogeneous bits/sections?  I was thinking, improve
statistics or similar.
But, honestly, at this point it doesn't really matter.  eta2 is so poorly
constrained anyways.

Simulation -- check model/fitting sensitivity
---------------------------------------------

Set up IPython notebook and letting it run now, w/ full model errors.
May be worth exploring large eta2 range, because the fit variation is MUCH more
pronounced at high eta2!  Consistent with what we are seeing in our results.
Running overnight to get full errors.

Paper writing
-------------

I am really hitting the limits of what I know (which is, not much).
I think, tomorrow I should spend a day just shredding through papers (and set
up full model fit errors to run, meanwhile).
Write up some material / take some notes.  Bring a notebook or buy a notebook,
whatever.


Tuesday 2014 September 23
=========================

Fitting to simulated data finished.
(15 hours to get errors for 6 regions! incl. parallelization.  geez)
Spent some time formatting data to send.

Reading/reviewing literature: goal is to review the relation between the
injected electron energy spectrum and the observed synchrotron spectrum.
Relations between diffusive/advective/acceleration time/length scales.
Most useful papers -- Reynolds 2008, Vink 2012 reviews.
X-ray brightness requires synchrotron cutoff, hence electron spectrum cutoff.

...

lots of algebra later.

Expression for cut-off energy verified, 
reviewed Parizot in detail, and several other papers.


Wednesday 2014 September 24
===========================

Paper writing.  Followed up a lot of references in Sean's paper.
Text, AAS abstract sent to Rob/Brian.


Thursday 2014 September 25
==========================

Summary
-------
* Data catchup.  Pipelining material on Tycho, Kepler, Cas A


Data catch-up: Tycho, Kepler, & Cas A
-------------------------------------

### Tycho

Double checked regions-all.reg, for other possible additions.
It seems like any possible options have already been ruled out, from 
previous iterations over `regions-2` or so.
The eastern limb, in particular, has spectral lines (see regions-2) even if I
only select the little rim section.
Southern regions are messed up by both low SNR and multiple filaments.

So, keep using regions-5.  Using FWHMs without cap, with background
subtraction, I generate specextract spectra for "up" regions.  Done around
noon (took 4 hrs).

Ran simple model fits + errors.
Now, running full model fits + full errors...

### Kepler

Going through `regions-1` -- inspect profiles, spectra, and spectrum fits.
Select/divvy up new regions, yielding 11 regions.  Kepler's ear was divided
into 7 regions, to get roughly equal counts in 2-7 keV band.

Generate 2-4, 4-7 keV counts images.

### Cas A

Picked a few regions to play with (get an idea of number of counts);
this is `regions-0`.  Generated some profile data + plots of profiles, to take
a look at on Monday -- no need to go further now.

Generate mosaic + counts images for 2-3, 3-4.5, 4.5-7, 7-9 keV.
(count images take a while for Cas A!)


Friday 2014 September 26
========================

Summary
-------
* Continued data catch-up.  Generated most new things needed.
* Manuscript updated / most recent Tycho plots/numbers (some still pending)

(mainly, waiting on full model errors to finish computing)

Data catch-up (cont): Tycho regions-5, Kepler regions-2
-----------------------------------------------------------

Updated `spec_linkbg.py` to work with DS9 fk5 files, for ease.
Updated `ds9projsplitter.py` to make region files w/ boxes directly.
(no need for separate script, since this is the only place where
`ds9proj2box.py` needs to be invoked)

### Kepler
Linked regions-1 spectra w/ backgrounds (forgot to do earlier).  Updated fits.
Reviewed my notes + new fits, it doesn't change our conclusions and region
picks.  Removing background does help w/ the very faint filaments.

For regions-2: Brian recommends throwing out regions where the profiles look
so bad/unusable.  Expectation is that physics is same everywhere...
(go ahead and do FWHMs, profile fits), especially since we already have
multiple samples of the ear filament.

Discard: region numbers 1, 2, 8, 10, 11.  Region number 9 is bordering on the
discard "threshold"; its highest energy band is a bit of a mess and doesn't fit
well.  The current black list is:

    blacklist = {1:['0.7-1kev'], 2:['0.7-1kev'],
                 8:['0.7-1kev'],
                 9:['4-7kev'],
                 10:['0.7-1kev', '4-7kev'],
                 11:['0.7-1kev']}

Generated all FWHMS (with/without cap, background subtraction) for model
fitting.

### Cas A
Brian notes that Fe K lines might be a problem.
I selected a bunch more filaments.  One big thing to consider: when using log
scaling, bright filaments get washed out; linear scaling does the reverse.
So don't forget to double check multiple image scalings/setups.

### Tycho
Tycho regions-5.  I linked upstream spectra to old bkg-2 spectra, which were
generated incorrectly.  Not perfect -- but it will do for now.
Just a minor issue but to be addressed.



Tycho paper
-----------

Removed old tables and figures (go to previous git commit to recover).
Most tables/figures now should be for Tycho regions-5 data.


(Week 18) Monday 2014 September 29
==================================

Summary
-------
* Meeting -- Tycho results, discussion, tables, paper fixes;
  Kepler simple model fits and regions; Cas A region selections
* Sent final draft AAS abstract to all
* Reviewed/added Tycho, Kepler SNR numbers (radio spectral index for Tycho)
* Kepler regions-2 simple model fits (first pass)
* Cas A regions profiles and FWHMs generated (from last Friday's new regions-0
  selections)

(set up specextract on bkg-2, to use all ObsIDs for consistency)

Monday meeting notes
--------------------
Also merged into (transient) to-dos

Tycho (paper disc. first priority): aim to send new version before Monday.
        lots of things to think about!  go through paper with Brian's comments
        Brian addressed a few questions / incorrect statements in the paper
        Make list of numbers/etc to update, fix, etc...
        Magnetic damping model is up next in line, very soon
        check error calculations on m_E etc...

Kepler: use new regions, run spectra, get full model numbers...
        otherwise looks quite fine

Cas A: take a look. select subset of regions, pull out NEARBY spectra.
       big issue is that cas A is so freaking bright that we have thermal
       contamination because of scattered light.  Geez.
       So we have to subtract or avoid that somehow.
       Data in 4-6.5 keV, 7+ keV would be ideal to side-step issue.
       Using photons below 4 keV -- gotta deal w/ contamination.


Kepler numbers
--------------

Looking up Kepler data values.  We need:
* distance to SNR (kpc)
* shock radius (arcsec)
* spectral (photon?) index s = 2 * alpha + 1
* shock velocity (cm/s)

I've gone ahead and updated SNR catalog w/ all relevant values.
They are NOT vetted with Rob and Brian, yet.

### Remnant distance

* Katsuda et al. (2008) give X-ray expansion measurements (6 yr baseline) for
  several regions around Kepler's SNR.  Regions 4/5 (SE filament) have
  expansions of 0.191, 0.206 arcsec/yr (+/- 0.03, at most, each).  Average:
  0.199 arcsec/yr = 3.77e8 cm/s +/- 0.57e8 cm/s velocity (stat+sys error is
  overestimate), assuming distance 4 kpc (4.71 +/- 0.71, if d = 5 kpc).

  3.3 kpc favored by Katsuda et al. (2008), matched Balmer filament
  measurements (Blair et al. 1991, Sankrit et al. 2005) to X-ray kinematics

* Vink (2008) gives 4200 km/s (assuming d = 4 kpc) for the bright SE filament,
  and states that "[this] is twice as fast as the shock velocities inferred
  from optical spectral and proper-motion studies in the [NW] region (Blair et
  al.  1991; Sankrit et al. 2005)..."; he favors a larger distance as being
  more consistent with the expected energetics/expansion of a Type Ia SN.

* Reynoso and Goss (1999) give lower bound of 4.8 +/- 1.4 kpc and upper bound
  of 6.4 kpc, by analysis of some absorption features -- I think the implicit
  assumption is that the HI velocities are indicative of location in/relative
  to the galactic plane?  Consider HI features behind/interacting with/in front
  of remnant.

* Patnaude et al. (2012) suggest ~5-6 kpc if expanding into ISM, and require
  >~7 kpc if expanding into a wind-blown bubble -- looking at energetics +
  hydrodynamic modeling.  A pure wind model is not enough; needs a small
  cavity to make the ionization ages work.
  (similar conclusion by Chiotellis et al. (2012))

  The TeV gamma ray non-detection implies distance > 6.4 kpc?  Quite extreme.
  (Aharonian et al. 2008) (what models do they consider to get this number?)

  A smaller distance (< 5-6 kpc) also favors a "subenergetic" SN.
  All of this is somewhat curious, that the different measurements and
  constraints are so discrepant.


Temporary conclusion -- take d = 5 kpc and run this by Rob/Brian...

### Shock radius

Shock radius = 1.5 arcmin (90 arcsec) from Green's catalog is fine.
In fact this is a little smaller than what I measure in DS9 -- especially since
we're interested in the protruding ear.  Chiotellis et al. (2012) give
radius of ~1.78 arcmin.  I measure radial dist to ears is about 2 arcmin,
radial dist to front of shock is about 1.8 arcmin.

Let's compromise and say ~1.9 arcmin.

### Shock velocity

For shock velocity -- taking the X-ray measured kinematics of Katsuda et al.
(2008) and scaling to 5 kpc gives 4.71e8 cm/s  (Vink (2008) would give 5.25e8
with d=5kpc).  The other 2 regions I have selected in regions-2 (10, 11)
are not sampled by Katsuda et al. (2008).  Reg-8 is closest to my region 10,
with an even faster expansion velocity, in fact; Reg-13 is closest to my region
11 but samples the slower northern rim, whereas region 11 falls on the ear of
Kepler (which is likely to be expanding more rapidly).

So I'll go ahead and assume shock velocity 4.71e8


### Radio spectral index (and then e- spectral/photon index)

The spectral index map (radio, 6-20 cm) of DeLaney et al. (2002) suggests that
at the nonthermal filaments, the radio spectral index is around -0.65 (but, the
N edge of the SE ear seems to be a bit steeper).  I will accept -0.64 from
Green's catalog for now.


Cas A region picks, FWHMs
-------------------------

Re-ran the radial profile generation scripts (forgot to keep the region files
ordered!  Messes everything up).  Ran FWHM fitting notebook and got decent
numbers, though in a number of cases the fitting code fails/crashes for
whatever reasons, not investigated.


Tycho radio spectral index
--------------------------

Working in reverse chronological order through various radio survey papers:
Sun et al. (2011), Sino-German 6cm Urumqi survey favor alpha = -0.58 +/- 0.02
(value in Green's catalog) over -0.65 from Kothes et al. (2006).  They fit many
radio flux measurements, from references:

  Langston et al. (2000)
  Fuerst et al. (1990b)  (ue is u-umlaut)
  Fanti et al. (1974)
  Bennett (1962)
  Conway et al. (1965)
  Kellermann et al. (1968, 1969)
  Pauliny-Toth et al. (1966)
  Kundu and Velusamy (1972)
  Green et al. (1975)
  Becker et al. (1991)
  Hurley-Walker et al. (2009)
  White and Becker (1992)
  Kothes et al. (2006)
  Horton et al. (1969)
  Reich et al. (1997)
  Bietenholz et al. (2001)
  Green (1986)

...

Okay I'm not even gonna challenge that.  The fit looks quite nice, though I
wonder what caused the jump from -0.65 of Kothes et al. to -0.58 (and, how
reliable the older data are).

* Kothes et al. (2006) actually got -0.61 from their Canada Galactic Plane
  Survey data alone, but pulled out Klein et al. (1979) to get -0.65 (funny,
  that doesn't seem to be in the seemingly exhausitve list of Sun et al.
  (2011))
* I'm not sure how Hurley-Walker et al. (2009) computed alpha in Table 6, but
  the fit plot for Tycho (Figure 2) gives alpha=0.58, which looks good.
* Katz-Stone et al. (2000) look for spectral index variation; find range
  -0.44 to -0.63, average -0.52



Tuesday 2014 September 30
=========================

Summary
-------
* Tycho paper updates/cleaning
* Attempt to submit AAS abstract foiled

Tycho data/paper updates / review
---------------------------------

### Data update

Error calculation finished last night around 9pm (so, mid-day Thurs to ~9p Mon
is about 4.5 days, 105 hours... geez).
Also computed full model fits with stderr only as backup.

Finished specextract / mergespectra on Tycho backgrounds-2.
Relinked spectra for regions-5 up.
Updated all tables (manually culling regions 21, 22; other regions untouched)
Updated all figures (SNR image, spectra, profile plots)

### Manuscript review

Propagating error through for `m_E`, turns out I was doing it wrong.  Woops.
Error should be:
    \[
        \delta m_E = \frac{1}{\ln(w_2/w_1)}
                     \sqrt{ \left(\frac{\delta w_2}{w_2}\right)^2 + 
                            \left(\frac{\delta w_1}{w_1}\right)^2 }
    \]
So even if `m_E` is zero, its error is certainly not!  Fixed in new tables.

Reviewed all of Brian's comments + added my own remarks.

Cleaned up and reorganized introduction/transport/diffusion text, cutting about
1 column of text.


Wednesday 2014 October 1
========================

Summary
-------
* Submitted AAS abstract
* Manuscript results/discussion, sketch and flesh out
* Updated manuscript spectra table/plots

Tycho manuscript
----------------

Reviewed and cleaned procedure/methods.  Cleaned results, outline and flesh out
some discussion.  Hit all comments, or at least integrated into text to be
fleshed out.

Quick chat with Brian about 1. where to place text on remnant distance
dependence (just 1-2 sentences...), and 2. radio rims.  Jack and Brian
have fresh VLA data (not yet run through VLA pipeline), but that likely won't
be tackled in earnest until December-ish, after new Chandra data comes in!

Remark: if we want to adapt the rim width model for radio filaments, we will
certainly have to delve back into the electron distribution calculation details.

Tycho data
----------

Run specextract on Tycho regions-5 "down" spectra, for completeness (to show
correct numbers and figures in draft -- I expect only FWHMs and fits to change
with regions-6).  Filled in tables/plots with correct downstream spectra (with
horrible lines and all).

Run `reproject_obs` on Tycho 750 ks data, to sanity check my mosaics and
generate new energy band images (since I don't have other programs/data to run
at the moment that are terribly pressing).  Should get same `merged_evt.fits`
that Brian sent, actually, I forgot I had that.  Woops.  Oh well -- will have
files on hand for consistency / meddling.

Just for fun, downloaded Hughes' 150ks Chandra observation of Tycho, from 2003.
I am wondering how much the rims may vary in time...
Attempted to run `chandra_repro` as usual, but gunzip fails with CRC check
failed on `acisf03837_000N003_evt1.fits.gz`.  Very strange.  This persists
after redownloading the data.


Thursday 2014 October 2
=======================

Magnetic damping code

Friday 2014 October 3
=====================

Filled out discussion and more conclusion (moving text around elsewhere to match)
Sent manuscript to Rob/Brian

Exploring effect of remnant distance on best fits
-------------------------------------------------

### Simple model

I compare simple model fits at d = 4 kpc instead of d = 3 kpc.
Looking at best fits to regions-5, FWHMS with subbkg (no cap).

B0 is basically unchanged.  eta2 scales as (4 kpc / 3 kpc)^2 = 1.77777; some
variation in the tables/numbers appears due to roundoff error.

At some very extreme values, B0 does vary slightly.  For example:

    |-----|----------------------|----------------------|
    |     | Region 16 (d = 4kpc) | Region 16 (d = 3kpc) |
    |     |----------------------|----------------------|
    | mu  | eta2        B0       | eta2        B0       |
    |-----|----------------------|----------------------|
    | 0   | 18705.6     8875.0   | 10362.2     8830.2   |
    | 1/3 | 10.4        860.2    | 5.9         860.2    |
    | 1/2 | 5.7         744.6    | 3.2         744.6    |
    | 1   | 2.2         618.2    | 1.2         618.2    |
    |-----|----------------------|----------------------|

At mu = 0, the eta2 ratio is 1.805, deviating from expected 1.777...
B0 ratio is 1.005.

But there are only a few cases of this, where B0 values actually change.
(in many cases though, the fitting routine hits the limit at eta2 = 10^5
or B0 = 10 mG)

### Full model

Basically the same, but the scaling is not quite so spot-on (as we might have
expected).  I have added my estimates of the deviation to the manuscript.
B0 values vary by about 1% or so; eta2 values vary by about 1-5% from the
expected 1.777... ratio.


Saturday 2014 October 4
=======================

Magnetic damping: planning out the chain of code modifications.

First get underlying code working (`models.width_cont` and everything below).
Then, deal with `models*.py` infrastructure to set-up damping calculations

    models_disp.py
    models_exec.py
    models.py

        width_cont **kwargs: idamp=False, damp_ab=0.05, damp_bmin=5.0e-6
                             pass to fefflen (only ever called by width_cont)

    FullEfflength_port.py

        I let emisx, distr (in FullEfflength_mod.f) do computations from
        B0/ab/Bmin separately -- no need to carry around B field array.
        e- distribution calculation uses z(x) = ... ugly thing

        fefflen args: idamp, ab, Bmin (ab, Bmin only used if idamp=True)
        emisx args: idamp, ab, Bmin

        pass idamp_flag, ab, Bmin to fullmodel.distr (idamp is an int)

    FullEfflength_mod.f

        distr(..., idamp, ab, Bmin)
        distrmlt1(..., idamp, ab, Bmin)
        distrmgt1(..., idamp, ab, Bmin)
        distrpohl(..., idamp, ab, Bmin)

        Also modified:
        Fullefflengthsub(..., idamp, ab, Bmin)

Finished updating code.  It looks like it works.  With Sean's default input I
get larger FWHMs but they definitely still look energy dependent.  Though, the
scale width was about 0.05 percent of remnant radius though.

I'm not sure whether to expect larger or smaller rims with damping.
Electrons radiate energy less effectively, so they get farther downstream
But the peak emission energy shifts down too, at smaller B
so they wouldn't radiate at that energy as much...

Oh well, it looks to work, will check it and all tomorrow.
(also, all the e- distr integrals' resolutions have to be checked)

Sunday 2014 October 5
=====================


B-damping code check
--------------------

Check my code against Sean's original B-damping code, ensure that we get the
same/similar numbers (for rim FWHMs) out

Phys parameters:
    eta = 1, mu = 1, Bmin = 5e-6, icut = 1
SN 1006 parameters:
    compratio = 4d0, v0 = 5d8/compratio, rs = 2.96e19,
    alpha = 0.6d0, s = 2*alpha + 1
Grid parameters:
    irmax, iradmax, ixmax = 1000, 400, 500

I generally take rminarc = 72, to match Sean's choice of rmin=0.92 (I assume
92% of 900 arcsec?).

### Sean's code

Compiled with `-O3` flag.  `rmin = 0.92`

* B0 = 50d-6, ab = 0.005
    1.00  keV   11.81  arcseconds   NaN  mnu
    2.00  keV   10.58  arcseconds  -.16  mnu
    4.00  keV    9.58  arcseconds  -.14  mnu
    8.00  keV    8.78  arcseconds  -.12  mnu
* B0 = 150d-6, ab = 0.005
    1.00  keV    7.85  arcseconds   NaN  mnu
    2.00  keV    6.84  arcseconds  -.20  mnu
    4.00  keV    5.98  arcseconds  -.19  mnu
    8.00  keV    5.26  arcseconds  -.19  mnu

* B0 = 50d-6, ab = 0.05
    1.00  keV   66.10  arcseconds   NaN  mnu
    2.00  keV   53.86  arcseconds  -.30  mnu
    4.00  keV   44.64  arcseconds  -.27  mnu
    8.00  keV   37.58  arcseconds  -.25  mnu
* B0 = 150d-6, ab = 0.05
    1.00  keV   15.84  arcseconds   NaN  mnu
    2.00  keV   12.10  arcseconds  -.39  mnu
    4.00  keV    9.65  arcseconds  -.33  mnu
    8.00  keV    7.92  arcseconds  -.28  mnu

### My modified code

Running code from `FullEfflength_mod.f` (NOT Python wrapper).
(tried Python wrapper code, but results failed too often. super sensitive to
rminarc value?...)

* B0 = 50d-6, ab = 0.005 (rminarc = 72 arcsec)
   1.00  keV:    16.92
   2.00  keV:    14.40
   4.00  keV:    12.42
   8.00  keV:    11.16
* B0 = 150d-6, ab = 0.005 (rminarc = 72 arcsec)
   1.00  keV:     9.72
   2.00  keV:     8.28
   4.00  keV:     7.02
   8.00  keV:     6.12

* B0 = 50d-6, ab = 0.05 (rminarc = 100 arcsec)
   1.00  keV:    78.50
   2.00  keV:    62.50
   4.00  keV:    50.25
   8.00  keV:    41.75
* B0 = 150d-6, ab = 0.05 (rminarc = 72 arcsec)
   1.00  keV:    16.56
   2.00  keV:    12.42
   4.00  keV:     9.72
   8.00  keV:     8.10

Something is up -- numbers are in the right ballpark, but rather off.  Too much
to be explained by resolution or similar (and, resolutions should be the same).

### Try running without energy cut-off (icut=0).  B0 = 150d-6, ab = 0.005
Sean's code:
    1.00  keV    9.43  arcseconds   NaN  mnu
    2.00  keV    8.71  arcseconds  -.11  mnu
    4.00  keV    8.06  arcseconds  -.11  mnu
    8.00  keV    7.63  arcseconds  -.08  mnu
My code:
    1.00  keV:   12.42
    2.00  keV:   11.34
    4.00  keV:   10.62
    8.00  keV:    9.90
Nope.  Besides, the expressions for cut-off energy look correct!

### Check resolutions?...

Woops, I used wrong resolutions in my Fortran set-up, fix that...
Sean's numbers, B0 = 150d-6, ab = 0.005
    1.00  keV    7.85  arcseconds   NaN  mnu
    2.00  keV    6.84  arcseconds  -.20  mnu
    4.00  keV    5.98  arcseconds  -.19  mnu
    8.00  keV    5.26  arcseconds  -.19  mnu
Now I get (B0 = 150d-6, ab = 0.005, cut-off enabled, rminarc=72):
   1.00  keV:    9.792
   2.00  keV:    8.424
   4.00  keV:    7.200
   8.00  keV:    6.264


### Disable magnetic damping

Set ab = 5d3, to trigger the flag and set `z(x) = x` (paper notation) in the e-
distribution calculation.  Synchrotron emissivity is still computed w/ damped
field.  But with such a large lengthscale it should not matter.

My code (rminarc = 72):
    1.00  keV:    15.768
    2.00  keV:    12.096
    4.00  keV:     9.648
    8.00  keV:     7.920
Sean's code:
    1.00  keV   15.55  arcseconds   NaN  mnu
    2.00  keV   11.95  arcseconds  -.38  mnu
    4.00  keV    9.58  arcseconds  -.32  mnu
    8.00  keV    7.85  arcseconds  -.29  mnu

Now, THAT might be close enough to chalk up to resolution error!
Try shrinking the resolution now.
My (Fortran) code w/ rminarc=20 gives:
    1.00  keV:    15.72
    2.00  keV:    12.08
    4.00  keV:     9.62
    8.00  keV:     7.96
Sean's code with rmin = 0.98 (~18 arcsec) gives:
    1.00  keV   15.64  arcseconds   NaN  mnu
    2.00  keV   12.01  arcseconds  -.38  mnu
    4.00  keV    9.58  arcseconds  -.33  mnu
    8.00  keV    7.90  arcseconds  -.28  mnu

Not amazing agreement.  We are using irmax = 1000, so resolution error is
about 20 arcsec / 1000 ~ 0.02 arcsec.  Moreover, my Fortran code implementation
uses the same stuff as Sean's code -- same FWHM calculation routine, etc.

AH, I just remembered, Sean is using splines to interpolate his intensity
profiles or whatever.  Try running my Python code w/ the FWHM finding routine,
to see if I get results closer to those of Sean's.

My Python code w/ rminarc=20, ab=5e3:
    1.00 keV: 15.77825
    2.00 keV: 12.12329
    4.00 keV: 9.67510
    8.00 keV: 8.00394

NOPE.  It got even farther away!  Why? -- ahh right, I remember.  Because the
FWHM calculation routine in Fortran code always takes the smallest FWHM value,
measured from grid points closest to peak / inside.

Well, at least FWHMs agree to order 0.1 arcsec (with FWHM values 8--16 arcsec)
I don't want to get into the nitty-gritty of the spline/integral calculations.
So let's say that this is acceptable.  Now, why do results with magnetic
damping enabled differ more strongly??

### Enable magnetic damping, use tight rmin/rminarc

Try again with ab = 0.005 (B0 = 150d-6, same resolutions, icut=1, etc)

Sean's code, rmin = 0.98 (basically the same, honestly)
    1.00  keV    7.87  arcseconds   NaN  mnu
    2.00  keV    6.88  arcseconds  -.19  mnu
    4.00  keV    6.03  arcseconds  -.19  mnu
    8.00  keV    5.33  arcseconds  -.18  mnu

My (Fortran) code, rminarc = 20:
    1.00  keV:    9.78
    2.00  keV:    8.38
    4.00  keV:    7.22
    8.00  keV:    6.28

Right, okay, we're back to this being off by 1--2 arcsec deal.
As a quick sanity check -- we are not invoking the advective solution so no
need to debug there (they are identical anyways).

But I look over the e- distribution calculations, and they are friggin
identical.  Seriously.  Well, let's take a look to be sure

### Check that calculated e- distributions agree

B0=150d-6, ab = 0.005, icut=1, irmax/iradmax/ixmax = 1000/400/500
Outputs to disttab2.dat, disttabmine.dat

Conclusion -- NUMBERS ARE IDENTICAL.  Okay.  So something is amiss.
Next step in line is to consider emissivities?

### Fix emissivity calculation in orig Fortran code

Realized that NOT been updating B field in this calculation (in Fortran code).
BUT, I am accounting for this in Python code and the Python code numbers are
even FARTHER off than Fortran...  Well, let's fix it and see:
(rminarc = 18)

   1.00  keV:    7.920
   2.00  keV:    6.912
   4.00  keV:    6.066
   8.00  keV:    5.364

Hallelujah.  We're now within about 0.05 of Sean's code.  Still not amazing,
but not bad.  I would expect agreement to +/- 0.018 if Sean did not use
splines.  I don't really want to go back and redo Sean's work without the
splines.  Now let's see why the Python code is so far off the mark.

CAUGHT: off by a sign.  Should be `exp(-1*(1-r)/ab)`.  Try the Python code!
(all same settings: rminarc = 18, ab = 0.005, Bmin = 5e-6, etc)
    1.00 keV: 7.96426
    2.00 keV: 6.96749
    4.00 keV: 6.11790
    8.00 keV: 5.40123

Better than before, but now this disagrees w/ Sean's code by abt 0.1, when we
expect them to agree within resolution (roughly) as well.
What about agreement with my own Fortran code?  Looks like generally offset by
~0.05 arcsec, seems promising.

### Retry all calculations with rminarc = 9 / rmin = 0.99

Make sure that the errors drop appropriately, if it is indeed tied to
resolution / spline usage / etc.

rminarc=9, rmin = 0.99; otherwise same settings as before.
ab = 0.005, Bmin = 5e-6, icut = 1, irmax/iradmax/ixmax = 1000/400/500

Expected resolution w/ irmax = 1000 is 9/1000 ~ 0.01 arcsec

* Aaron's Python code
    1.00 keV: 7.96428
    2.00 keV: 6.96751
    4.00 keV: 6.11794
    8.00 keV: 5.40130
* Aaron's Fortran code (modified version of Sean's original)
    1.00  keV:    7.920
    2.00  keV:    6.930
    4.00  keV:    6.084
    8.00  keV:    5.373
* Sean's Fortran code (for B-damping)
    1.00  keV    7.87  arcseconds  -.12  mnu
    2.00  keV    6.89  arcseconds  -.19  mnu
    4.00  keV    6.05  arcseconds  -.19  mnu
    8.00  keV    5.33  arcseconds  -.18  mnu

But, the errors are holding constant at about 0.04 arcsec, in stepped fashion.

Try one more calculation, all same except let ab = 0.5 instead.  Much weaker
damping -- magnetic field drops by ~2% over 9 arcseconds; about 10% over 45
arcseconds.  Let rminarc = 36, rmin = 0.96.  Then...

* Aaron's Python
    1.00 keV:   15.83766
    2.00 keV:   12.15671
    4.00 keV:    9.69488
    8.00 keV:    8.01593
* Aaron's Fortran
    1.00  keV:  15.768
    2.00  keV:  12.096
    4.00  keV:   9.648
    8.00  keV:   7.956
* Sean's Fortran
    1.00  keV   15.70  arcseconds   NaN  mnu
    2.00  keV   12.02  arcseconds  -.38  mnu
    4.00  keV    9.58  arcseconds  -.33  mnu
    8.00  keV    7.92  arcseconds  -.27  mnu

Mismatch between Sean's Fortran and my Fortran is 0.03 to 0.07 arcsec
Mismatch bewteen Sean's Fortran and my Python  is 0.1 to 0.13 arcsec

So, again, we see this stepped error.  Resolution error here should be about
0.04 arcsec, compared to before.

The error here, then, would be of order 1%.  But this is not great, because my
resolution error bounds are at *most* 1% (and typically 0.1%).  So this would
be a dominant error term.
(and, I don't really have a way to validate whether my code is right, or
Sean's, or whatever)




(Week 19) Monday 2014 October 6
===============================

**EDIT 2014 October 30: WARNING, equations/numbers have errors!
(missing factor of 100e-6 G, code given has wrong exponent hidden inside).
Do NOT use provided values of eta2**

Summary
-------
* Discussion w/ Rob (misc. small things, travel, synchrotron roll-off)
* Looked over paper, added small comments/fixes
* Reviewed roll-off freq arg, derived eta2 prediction.  Tested w/ srcut

Mailed in AAS junior membership form w/ dues.

Meeting with Rob
----------------
* Magnetic damping -- don't spend too much time on it.  Just run it through...
* On using synchrotron roll-off -- try it with srcut.
  Kind of hard because we don't have that kind of spatially resolved radio
  spectroscopy to constrain radio flux, but give it a try... (probably won't
  pan out? kind of tangential)
* Roll-off -- Rob pulled up Tanaka et al. (2008, ApJ), Suzaku paper on
  RX J1713.7-3946 spectrum in hard & soft X-rays
* Small theory details (why Bohm limit?  why D propto E^2? etc), check w/Steve
* Email Dave Holdridge about travel (done)


Synchrotron roll-off and diffusion
----------------------------------

Idea: the synchrotron roll-off could give an independent constraint on
diffusion in each of our regions/filaments.  DSA gives a prediction for the
synchrotron roll-off in terms of diffusion coefficient only.  Whereas we
estimated magnetic field and diffusion from the rim widths, set by diffusion
affecting electrons propagating downstream of the shock, this estimates
diffusion strictly from the DSA process at/near the shock.  If diffusion
coefficients vary behind the shock these values could be quite different (would
make sense in context of a magnetic damping model too).

### Derivation
(typing up paper notes)

Assume that a DSA model explains the observed e-/synchrotron spectra
cut-off/roll-off.  We estimate e- cut-off energy by equating
`$\tau_{\mt{accel}} = \tau_{\mt{synch}}$`.  This gives (Parizot et al., 2006):

    \frac{3r}{r-1} \frac{r D_d + D_u}{v_s^2} = \frac{1}{b B^2 \Ecut}

(\emph{I don't know if this equality was designed to exactly obtain the electron
cutoff energy, or is just a rough/scaling relation -- pull up Drury (1983) to
check this})
Following through the arguments/assumptions in Parizot et al. (2006) yields
(using `$13.3 \unit{erg} = 8.3 \unit{TeV}$`):

    \Ecut =
        &\left( 13.3 \unit{erg} \right)^{2/(1+\mu)}
        \left( \frac{B_0}{100 \muG} \right)^{-1/(1+\mu)} \nonumber \\
        &\times
        \left( \frac{v_s}{10^8 \unit{cm\;s^{-1}}} \right)^{2/(1+\mu)}
        \eta^{-1 / (1+\mu)} .

Now, we can convert this cut-off e- energy to a roll-off synchrotron energy by
invoking the delta function assumption `$\nu_{\mt{cut}} = c_m \Ecut^2 B$` with
`$c_m = 1.822 \times 10^{18}$`.
(\emph{Again, I don't know if this is an equality or just a scaling.  We'd need
to relate the electron spectrum cut-off to the derived synchrotron spectrum
cut-off, as is done (I think) by Zirakashvili and Aharonian (2007).})

    \nu_{\mt{cut}} = c_m
        \left( 13.3 \unit{erg} \right)^{4/(1+\mu)}
        \left( \frac{B}{100\muG} \right)^{-(1-\mu)/(1+\mu)}
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{4/(1+\mu)}
        \eta^{-2/(1+\mu)}

In the case of `$\mu = 1$` this yields:

    \nu_{\mt{cut}} = (0.133 \unit{keV} / h) 
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{2}
        \eta^{-1}

Now, for `$\mu \neq 1$` this is no longer independent of magnetic field.  But,
there is hope.  Rewrite this in terms of `$\eta_2$` using
`$\eta = \eta_2 E_2^{1-\mu}$` and write:

    \nu_{\mt{cut}} = c_m
        \left( 13.3 \unit{erg} \right)^{4/(1+\mu)}
        \left( E_2 \right)^{-2(1-\mu)/(1+\mu)}
        \left( \frac{B}{100\muG} \right)^{-(1-\mu)/(1+\mu)}
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{4/(1+\mu)}
        \left( \eta_2 \right)^{-2/(1+\mu)}

With `$E_2 = \left( (2 \unit{keV}/h) / (c_m B) \right)^{1/2}$`, or
`$E_2 = \left( 0.2657 \unit{erg^2\;G} / B \right)^{1/2}$`, we obtain:

    \nu_{\mt{cut}} = c_m
        \left( 13.3 \unit{erg} \right)^{4/(1+\mu)}
        \left( 0.2657 \unit{erg^2\;G} / B \right)^{-(1-\mu)/(1+\mu)}
        \left( \frac{B}{100\muG} \right)^{-(1-\mu)/(1+\mu)}
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{4/(1+\mu)}
        \left( \eta_2 \right)^{-2/(1+\mu)}

And, the magnetic field terms cancel!  Hey presto.

### Result (synchrotron roll-off and eta2)

Our final usable result is:

    \nu_{\mt{cut}} = c_m
        \left( 13.3 \unit{erg} \right)^{4/(1+\mu)}
        \left( 0.2657 \unit{erg^2\;G} \right)^{-(1-\mu)/(1+\mu)}
        \left( 100 \muG \right)^{(1-\mu)/(1+\mu)}
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{4/(1+\mu)}
        \left( \eta_2 \right) ^{-2/(1+\mu)}

Now we have a prescription for estimating `$\eta_2$` from the
synchrotron cut-off frequency and the shock velocity, no other info needed.
You can verify that this reduces correctly when mu=1.
Use this Python function for one-off interactive calculations:

    from __future__ import division
    
    def cut_eta2(nu_cut, vs_sc, mu):
        """Compute eta2 from synchrotron roll-off frequency, assuming that
        e- roll-off occurs when synchrotron cooling timescale (~1/(B^2 E))
        balances DSA timescale (~D/v^2)

        Input:
            nu_cut = synchrotron roll-off in Hz (from srcut or similar)
            vs_sc = shock velocity scaled by 10^8 cm/s
            mu = diffusion-energy scaling exponent (mu=1 for Bohm)
        Output:
            eta2 estimate
        """

        nu_pred = ( 1.82e18 * (13.3)**(4./(1+mu)) * 0.2657**(-(1.-mu)/(1.+mu))
                    * (1e-4)**(2./(1.+mu)) * (vs_sc)**(4./(1+mu)) )
        eta2 = (nu_pred / nu_cut)**((1.+mu)/2.)
        return eta2

If the cut-off frequency is less than 2 keV/h, eta2 will increase with mu.
If the cut-off frequency is greater than 2 keV/h, eta2 decreases with mu.

Why?  eta2 is independent of mu at 2 keV so it doesn't matter there.
IF the cut-off freq occurs below 2 keV, larger mu would DECREASE eta2; thus
eta2 must increase to reproduce the correct diffuion at the cut-off energy.
And vice versa -- if cut-off freq occurs above 2 keV, larger mu would INCREAse
eta2 and thus eta2 must drop to get the right diffusion.

This is specifically an artifact of our choice of fiducial energy at 2 keV.

Then, the more interesting question is -- why does eta2 decrease with mu
in our modeling efforts?  Is this at all meaningful?

Compute eta2 from srcut frequency for two regions
-------------------------------------------------

Using `phabs*srcut` model

### Region 6 from regions-5

    ========================================================================
    Model phabs<1>*srcut<2>  Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.646815     +/-  2.05401E-02  
       2    2   srcut      alpha               0.580000     frozen
       3    2   srcut      break      Hz       6.95027E+16  +/-  6.72192E+15  
       4    2   srcut      norm                4.02326E-02  +/-  4.76079E-03  
    ________________________________________________________________________
    
    
    Fit statistic : Chi-Squared =         213.44 using 234 PHA bins.
    
    Test statistic : Chi-Squared =         213.44 using 234 PHA bins.
     Reduced chi-squared =        0.92398 for    231 degrees of freedom 
     Null hypothesis probability =   7.903088e-01

The size of region 6 (up slice) is 0.0496 * 0.214 arcmin^2 = 0.0106 arcmin^2.
Compare Tycho's size is 4 arcmin radius, or 50.27 arcmin^2.  Ratio of areas is
0.00021; ratio of norm/total radio flux is (0.04 Jy)/(56 Jy) = 0.0007.
Does not seem unreasonable, since this is likely a brighter part of the
remnant.  So looks fine to me.

Assume shock velocity `3.9e8 cm/s * 3/2.3 = 5.1e8 cm/s` and accept srcut
roll-off freq 6.95e16 Hz (0.287 keV).  This gives:

    mu      eta2
    -------------
    0        4.57
    1/3      6.31
    1/2      7.42
    1       12.05
    1.5     19.57
    2       31.78

Verdict: vaguely plausible?  I suspect this would give more consistent magnetic
field estimates.


### Region 17 from regions-5

    ========================================================================
    Model phabs<1>*srcut<2>  Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.590624     +/-  2.39908E-02  
       2    2   srcut      alpha               0.580000     frozen
       3    2   srcut      break      Hz       7.29093E+16  +/-  8.86814E+15  
       4    2   srcut      norm                2.81425E-02  +/-  4.09850E-03  
    ________________________________________________________________________
    
       Using energies from responses.
    
    Fit statistic : Chi-Squared =         173.79 using 190 PHA bins.
    
    Test statistic : Chi-Squared =         173.79 using 190 PHA bins.
     Reduced chi-squared =        0.92937 for    187 degrees of freedom 
     Null hypothesis probability =   7.470014e-01

    XSPEC12>error 1,3,4
     Parameter   Confidence Range (2.706)
         1     0.554489     0.633984    (-0.0361356,0.0433598)
         3  5.91195e+16  8.80831e+16    (-1.37899e+16,1.51737e+16)
         4    0.0225397    0.0366237    (-0.00560276,0.0084812)

We take roll-off freq 7.29e16 Hz (0.302 keV) with
shock velocity 3300 cm/s x 3/2.3 = 4300 cm/s.  This gives:

    mu      eta2
    -------------
    0        3.17
    1/3      4.35
    1/2      5.09
    1        8.17
    1.5     13.10
    2       21.03

Again, verdict, vaguely plausible?  This procedure will yield very similar
numbers all around the remnant.   The only things changing are cut-off
frequency (varying by 5% maybe) and shock velocity (varying by factor of 2).


Tuesday 2014 October 7
======================

Summary
-------

* Ran Kepler regions-3 through pipeline
* Manually played with B-damping FWHMs, check weird behavior
* Set up calculation (for one set of data, mu=1) to run overnight w/ several
  scale lengths and a range of initial parameters (check how chisquared space
  looks?..)

Kepler regions-3
----------------

Short diversion.  Looking at Kepler notes from Sept 26: Brian wanted to
throwing out regions where profiles bad, even at lowest bands.  The blacklist
would be 1,2,8,9,10,11 -- but that seems a bit too many.  It would leave us
with just 5 regions (3-7) all on Kepler's ear.

Regions 9, 10 must go -- their bands with low counts look terrible.
Region 1 is iffy... to be safe I will throw it out.

So regions-3 we will have 2-7, 8, 11 (total 8 regions to get profiles/fits)

Set up profiles/FWHMs for Kepler as usual (same settings, Hanning window length
13 for smoothing), ran specextract.  Linked backgrounds and ran spectral fits,
same procedure as Tycho.  Ran simple model fits.  Ready to go on full model...


B-damping code
--------------

Filling in rest of framework for running fits...

    tycho-fit-tables.ipynb
    
        Iterate over data with models_exec.Fitter generator
        models_disp.build_dataf to generate fitting function
        models_disp.generate_fits to run and save fits

    models_disp.py

        build_dataf passes **fit_kws to Fitter.get_fit(...),
                           **err_kws to Fitter.get_errs(...)

        The rest of functions don't affect fitting/model calculations,
        just run or display fits (using function from build_dataf).
        build_dataf is mainly to extract and parse data from fits and error
        computations.  Fits/errors are run from one-liners

    models_exec.py

        Fitter.get_fit(...) passes **fit_kws to
            models.simple_fit, models.full_fit

            This method adds 'scale_covar':False to **fit_kws if not already
            given.  For full fit, if 'B0'/'eta2' not added, method adds best
            grid 'B0'/'eta2' from stored grid using self.grid_scan(mu).

    models.py

        simple_fit, full_fit accept kws:
            eta2, B0, mu_free, eta2_free, B0_free, **lmf_kws
        full_fit also accepts the kwarg 'model_kws'

            **

            (where **lmf_kws are settings passed to lmfit.minimize, and also
            scipy.leastsq)

As a small aside, I removed `Fitter.fitter_simp` and `Fitter.fitter_full` since
they subsetted the functionality of `Fitter.get_fit` but were less useful.

For manual fitting / getting numbers to play with, we have

    models_exec.py

        Fitter.width_full(..., **kwargs) passes kwargs to models.full_width,
        which passes **kwargs straight to models.width_cont.


Example B-damping calculation/numbers
-------------------------------------
With Bmin = 5e-6 (Sean's default).

Best fit with no damping is eta2 = 3.7 (+2.4/-1.3), B0 = 208 (+17/-12), with
chi-squared = 25.

    ab      eta2    B0      chisqr
    0.5     3.86    210     25.0
    0.05    2.96    195     25.5
    0.007   0.37    29.1
    0.005   0.03    33.9    25.79
    0.002   failed miserably (stalling at FWHMs too small)
    0.001   failed miserably

Examples of fits where eta2 is fixed

    ab      eta2    B0      chisqr
    0.005   0.01    35.6    -       (didn't record)
    0.005   1.00    22.92   85.6

At fixed eta2, I'm getting the impression that the widths are no longer nearly
as sensitive to B0 as before, and MUCH smaller B0s are required just to get
anywhere!  How should I best get a handle on parameter space?


Manual B-damping calculations (log of testing)
----------------------------------------------

Testing some numbers... determined that decreasing eta2 can increase rim widths
in the damped model, due to two effects:
* Cut-off energy increases with smaller eta2, allowing larger rims (more high
  energy e- available)
* Without Ecut, rims still widen w/ smaller eta2 (which was very perplexing).
  I think the answer lies in the diffusion equation:

    v df/fx - d/dx(D df/dx) + ... = v df/dx - D d/dx(df/dx) - dD/dx df/dx

  The last term, dD/dx df/dx, causes advection to slow down in the presence
  of a strong _diffusion coefficient_ gradient!  This is really trippy and I
  don't have an intuitive understanding of why this is.

(I discarded my log, but the same effect persisted in Sean's original Fortran
code -- so it's not just me)

A bigger issue now is that the FWHM calculation fails on some (large) values of
rminarc.  See the below code logging...

mu = 1, Tycho Region 1, irmax = 200 fixed.

First, just vary eta2 w/ other settings fixed (default settings)

    eta2    B0      rminarc icut    idamp   ab      FWHMs (1 keV, 2 keV, 3 keV, 4.5 keV)
    -------------------------------------------------------------------------------------------------
    0.03    34      10      1       1       0.004   [ 7.2240338   5.89908341  5.2959738   4.77683576]
    1       34      10      1       1       0.004   [ 3.29273975  2.91199758  2.73110823  2.57283302]
    10      34      10      1       1       0.004   [ 2.34055386  2.13230373  2.01590829  1.90183956]

Check if effect persists, when cut-off energy not introduced

    0.03    34      10      0       1       0.004   [ 240.        240.        8.45320352  6.82156553]
    1       34      10      0       1       0.004   [ 4.4062853   3.96766419  3.79779506  3.67446985]
    10      34      10      0       1       0.004   [ 3.58824847  3.52330823  3.50092131  3.48572739]

The code is also getting weird results for large rminarc, which I do NOT expect.

    0.03    34      10      1       1       0.004   [ 7.2240338   5.89908341  5.2959738   4.77683576]
    0.03    34      30      1       1       0.004   [ 240.        5.89989288  5.29563315  4.77711038]

The scourge is found when I disable adaptive rminarc calculation.  When
adaptive rminarc is on, it drops rminarc too small (overshooting), so that the
FWHM isn't found and 240. is reported.  (fix) denotes `irad_adapt=False`

    0.03    34      30(fix) 1       1       0.004   [-11.84010424 5.81935857  5.23547149  4.7287726 ]

Quick sanity check that all is normal, when damping is off

    1       100     60      1       0       -       [ 21.93518981  16.72890036  14.53880802  12.8070013 ]
    10      100     60      1       0       -       [ 33.24745697  28.39278335  26.15818132  24.23900033]
 

Next steps
----------

Set up some small B-damping calculations to run manually overnight.
I don't really know what are reasonable starting parameter values, yet.
So I'm entering my guesses in by hand... as we progress I'll remedy this.

Tomorrow: run things by Rob/Brian, hunt down rminarc bug again...
          tackle resolution error
          explore fitting outliers (in non-damping code)


Wednesday 2014 October 8
========================

Summary
-------
* Preliminary B-damping results (email, talk w/ Brian)
* Debug B-damping FWHMs
* B-damping theory -- diffusion gradient, eff. veloc, lengthscales from YPL '05
* Prep for B-damping tables
* Plan of attack: resolution errors, B-damping, downstream spectra,
    best fit mE values, srcut-based model fits


Running into Rob on a Wednesday morning
---------------------------------------

Rob's main observation on paper draft -- need to add fits that either excise or
fit Si/S lines, to show that it's just the lines and there's still a power law
underneath.  So do that.

Assess manual B-damping calculations
------------------------------------

First look at results -- goodness.  We can get pretty good fits with
moderate/strong damping.  This is just for Region 1 of Tycho.  Generated
results for Region 15 too.

Note Region 15 best fits favor somewhat weaker energy-dependence because error
on lowest energy band is so large, and it disproportionately does-not-influence
the fit.  In fact, this may be an issue for many of our fits...

Actual results are given at the bottom of this day's notes.


B-damping profile FWHM observations & debugging (cont.)
-------------------------------------------------------
(cont. from yesterday that is)

Fixed negative FWHM results -- code searched for crossings on grid of intensity
values.  Spurious crossings due to local minimum in intensity, just behind
peak, were misidentified as real and yielded negative FWHMs.

rminarc errors (box apparently too small) also arise when peak simply does not
have a FWHM!  Behind the peak there is too much emission, staying above the
peak half max.  Again, debugging with extreme cases, but this is an interesting
observation.

Weird behaviors observed with B-damping intensity so far, summarized:

* Stronger diffusion can cause rims to _narrow_ due to coefficient gradient
  (see yesterday's notes)
* Strong emission background behind peak can create a castle parapet-like
  morphology (which is, similar to what is seen in some locations -- but it's
  uncertain whether that's due to thermal emission)
* Associated w/ strong emission background, sometimes, is this local minimum
  behind the peak.  I.e., emission increasing towards remnant center?!
  (or, possibly a second hump/peak behind the main peak, I don't know)

These could be interesting/useful, if fitting entire profiles (as others have
tried)!  But farther downstream of shock we're more likely we are to get hit with
other processes / going-ons.  Who knows what happens to the B-field and
diffusion coefficient, out there.

### More derivation -- effective advection velocity?

See some notes/derivation I scribbled down.  If this very simple description of
advection speed dropping due to diffusion coefficient gradient is correct,
there is a VERY strong variation in "effective advection velocity" at the
shock, persisting over several damping lengthscales.  My very quick
calculations even show this "effective velocity" going strongly negative, to
-1/2 x shock speed, before rising up again.

(this then would be consistent with this idea of a (very weak / barely present)
local minimum just behind rim profiles, or so)

A caveat is that our transport equation was for a plane shock and assumed
constant downstream veloc. etc... so there are almost certainly more
terms, effects, etc.



Discussion w/ Brian
-------------------

Our results/presentation may change dramatically.
We may not even need to present our fit tables if it turns out that we just
have this negative result, that we can't distinguish anything much.

B-damping tables/results first, to further vet this.
Brian asked something like -- what's a physically reasonable number?  What kind
of lengthscale should we consider?  See my notes below for answer (basically
the range of ab we're looking at is okay, for Tycho).

Brian will send srcutlog, an XSPEC "local model" -- better for fitting, not
included w/ XSPEC.  The break freq. / normalization might be somewhat
correlated if both left free, per Steve (in Brian's recollection).
The srcut result would be useful.  At this point the model is so iffy that it's
difficult to trust the results for anything but B0 values.

Results for meeting on Monday morning as usual (or, even earlier)

Brian can get to paper on Friday (XMM-Newton proposal deadline is something
like 7am EDT (noon UTC)).  Revisions to paper will have to wait on discussion
of results anyways.


B-damping lengthscales
----------------------

### Rough expectation

Pohl/Yan/Lazarian 2005 (PYL) suggest a typical lengthscale of 10^16, 10^17 cm
(considering 3 different types of energy cascades that could damp the field).

Tycho's radius is 240 arcsec = pi/2700 radians.  Assuming d = 3kpc, its
radius is then 1.077 x 10^19 cm, or 3.5 pc.  Then the (scaled) damping scale
length ab would be 0.001-0.01x shock radius for the range of numbers
10^16-10^17 cm stated by PYL.

(by comparison, SN 1006 radius is 900 arcmin w/ d = 2.2kpc, so radius
1.77x 10^20 cm would give ab ~ 0.00006 to 0.0006 (6e-5 to 6e-4)

### Compute with SNR specific parameters (predictions for Tycho)

The damping lengthscale also scales with shock velocity, sqrt(ISM density),
turbulent wavelength, and a measure of B fluctuation (variance or something).

For us, U ~ 5000 km/s.  n is ~0.1 to 0.5, but sqrt(n) will give something
closer to 1.  The turbulent wavelength lambda is uncertain, but let's assume
it is simply the e- Larmor radius (I have no idea why).

Synchrotron photons at 1-9 keV are radiated by 23-68 TeV e- in a 100 muG field.
In a 20 muG field we get 52-150 TeV e- energies.  From equation (4) of PYL, we
get Larmor radii of 0.7--2 x 10^15 cm for 100 muG field, 8--23 x 10^15 cm for a
20 muG field.

Finally, for B fluctuation (delta B), let's assume 100 microG.  PYL cite the
Bell-Lucek streaming instability for this, relevant for a young SNR -- but I'm
not sure precisely how young.

Then, e.g., for Tycho we have:

    (Kolmogorov cascade)
    l_d(B=100 muG) = (10^16 cm) * 5 * sqrt(0.1) * (0.7 to 2) = 1e16 to 3e16
    l_d(B=20 muG) =  (10^16 cm) * 5 * sqrt(0.1) * (8 to 23) = 1.3e17 to 4e17

    (YL fast-mode cascade, assuming L = shock radius)
    l_d ~ (10^16 cm) * sqrt(3.5/3) * sqrt(0.7 to 2) = 0.9e16 to 1.5e16
    l_d ~ (10^16 cm) * sqrt(3.5/3) * sqrt(8 to 23) = 3e16 to 5e16

    (Alfven mode cascade)
    Same as for YL, but multiplied by 5

Slightly more specific numbers for Tycho, assuming turbulent wavelength =
Larmor radius, assuming some velocity factors are unity, assuming L is SNR
radius, assuming typical field fluctuation of 100 muG.

* If B = 100 muG, damping length range is   1--8 x 10^16 cm
* If B =  20 muG, damping length range is 0.3--4 x 10^17 cm

*Conclusion*: ab of about 0.001 to 0.04 x shock radius is reasonable for Tycho.


Manual B-damping fit results
----------------------------

ab = damping lengthscale scaled to shock radius (ab=0.005 is 1.2 arcsec)
Ffit? indicates whether eta2 was free parameter or not
First set of eta2/B0 columns is (hand-input) initial guess for fit
Second set of columns contains best fit parameters
chisqr is chisqr, not reduced/scaled by DoF.
ID is an ID associated with plots / saved output

Output fits were saved to:

    full-std_err-damping-141007-test-22-fobj.pkl
    full-std_err-damping-141008-reg15-test-22-fobj.pkl

(the unlabeled one is Region 1, of course)

### Tycho Region 1

    ab      Ffit?   eta2    B0          eta2    B0      chisqr      ID
    ------------------------------------------------------------------
    0.500   True    3.700   208.000      3.854  209.675      25      1
    0.050   True    3.000   195.000      2.958  195.038    25.5      2
    0.040   True    3.000   210.000      2.308  183.488    25.7      3
    0.030   True    2.000   200.000      1.381  162.803    26.1      4
    0.020   True    2.000   200.000      0.375  128.534    27.6      5
    0.010   True    1.000   100.000      0.993   32.078    32.6      6
    0.009   True    1.000   100.000      1.278   27.554    28.1      7
    0.008   True    1.000   100.000      1.630   24.394    24.3      8
    0.007   True    0.500   100.000      0.002   44.256    67.6      9
    0.006   True    0.500   30.000       0.090   32.905    25.6     10
    0.005   True    0.030   34.000       0.032   33.914    25.8     11
    0.004   True    0.030   34.000       0.011   34.275    25.1     12
    ------------------------------------------------------------------
    0.500   False   1.000   210.000      1.000  180.958    35.5     13
    0.050   False   1.000   210.000      1.000  177.284    32.2     14
    0.040   False   1.000   210.000      1.000  172.488    29.7     15
    0.030   False   1.000   210.000      1.000  160.596    26.7     16
    0.020   False   1.000   210.000      1.000  121.947    31.2     17
    0.010   False   1.000   100.000      1.000   32.034    32.6     18
    0.009   False   1.000   100.000      1.000   28.812    28.2     19
    0.008   False   1.000   100.000      1.000   26.436    24.4     20
    0.007   False   1.000   100.000      1.000   24.733    24.8     21
    0.006   False   1.000   100.000      1.000   23.591    38.4     22
    0.005   False   1.000   23.000       1.000   22.918    85.6     23
    0.004   False   1.000   23.000       1.000   22.628     201     24

    0.500   False   0.100   200.000      0.100  168.458    71.1     25
    0.050   False   0.100   200.000      0.100  168.602    66.5     26
    0.040   False   0.100   200.000      0.100  165.932    60.4     27
    0.030   False   0.100   200.000      0.100  158.969    48.4     28
    0.020   False   0.100   200.000      0.100  136.747    30.6     29
    0.010   False   0.100   100.000      0.100   51.367    41.5     30
    0.009   False   0.100   50.000       0.100   44.098    38.8     31
    0.008   False   0.100   50.000       0.100   38.788    34.2     32
    0.007   False   0.100   50.000       0.100   35.036    28.5     33
    0.006   False   0.100   50.000       0.100   32.527    25.6     34
    0.005   False   0.100   27.000       0.100   31.040    39.5     35
    0.004   False   0.100   23.000       0.100   23.000   4e+06     36

    0.500   False   10.000  220.000     10.000  245.438      28     37
    0.050   False   10.000  220.000     10.000  230.338    29.5     38
    0.040   False   10.000  220.000     10.000  218.482    31.4     39
    0.030   False   10.000  200.000     10.000  189.565    36.2     40
    0.020   False   10.000  100.000     10.000   95.351    47.5     41
    0.010   False   10.000  100.000     10.000   21.390    29.5     42
    0.009   False   10.000  100.000     10.000   19.966    26.1     43
    0.008   False   10.000  100.000     10.000   18.874    24.9     44
    0.007   False   10.000  40.000      10.000   18.068    30.9     45
    0.006   False   10.000  40.000      10.000   17.521      55     46
    0.005   False   10.000  27.000      10.000   17.200     119     47
    0.004   False   10.000  23.000      10.000   17.080     254     48

### Tycho Region 15

    ab      Ffit?   eta2    B0          eta2    B0      chisqr      ID
    ------------------------------------------------------------------
    0.500   True    3.700   208.000     67.926  899.200    19.8      1
    0.050   True    3.000   195.000     89.259  940.063    20.2      2
    0.040   True    3.000   210.000     84.523  922.534    20.3      3
    0.030   True    2.000   200.000     81.244  902.027    20.6      4
    0.020   True    2.000   200.000     58.524  815.319    21.1      5
    0.010   True    1.000   100.000     28.823  601.958    20.7      6
    0.009   True    1.000   100.000     23.912  552.525    20.4      7
    0.008   True    1.000   100.000     15.754  485.146    20.2      8
    0.007   True    0.500   100.000      7.766  406.193      20      9
    0.006   True    0.500   30.000      69.771  108.082    20.3     10
    0.005   True    0.030   34.000      46.996   48.548    20.3     11
    0.004   True    0.030   34.000      23.054   33.049    20.2     12
    ------------------------------------------------------------------
    0.500   False   1.000   210.000      1.000  426.463    70.3     13
    0.050   False   1.000   210.000      1.000  429.314    72.1     14
    0.040   False   1.000   210.000      1.000  429.605    72.2     15
    0.030   False   1.000   210.000      1.000  429.574    71.9     16
    0.020   False   1.000   210.000      1.000  427.440    69.6     17
    0.010   False   1.000   100.000      1.000  402.638    51.3     18
    0.009   False   1.000   100.000      1.000  393.183    46.1     19
    0.008   False   1.000   100.000      1.000  379.205    39.7     20
    0.007   False   1.000   100.000      1.000  357.710    32.4     21
    0.006   False   1.000   100.000      1.000  322.264    25.1     22
    0.005   False   1.000   23.000       1.000  257.941    20.3     23
    0.004   False   1.000   23.000       1.000  133.117    21.1     24

    0.500   False   0.100   200.000      0.100  388.367     156     25
    0.050   False   0.100   200.000      0.100  392.284     161     26
    0.040   False   0.100   200.000      0.100  392.954     161     27
    0.030   False   0.100   200.000      0.100  393.722     162     28
    0.020   False   0.100   200.000      0.100  393.641     159     29
    0.010   False   0.100   100.000      0.100  380.135     130     30
    0.009   False   0.100   50.000       0.100  374.472     121     31
    0.008   False   0.100   50.000       0.100  366.088     107     32
    0.007   False   0.100   50.000       0.100  353.311    89.9     33
    0.006   False   0.100   50.000       0.100  332.992    67.8     34
    0.005   False   0.100   27.000       0.100  298.251    42.9     35
    0.004   False   0.100   23.000       0.100   23.000 7.91e+0     36

    0.500   False   10.000  220.000     10.000  597.601    22.7     37
    0.050   False   10.000  220.000     10.000  597.754    23.5     38
    0.040   False   10.000  220.000     10.000  596.889    23.7     39
    0.030   False   10.000  200.000     10.000  594.471    23.9     40
    0.020   False   10.000  100.000     10.000  585.725    23.8     41
    0.010   False   10.000  100.000     10.000  522.486    21.4     42
    0.009   False   10.000  100.000     10.000  499.541    20.8     43
    0.008   False   10.000  100.000     10.000  465.474    20.2     44
    0.007   False   10.000  40.000      10.000  411.710      20     45
    0.006   False   10.000  40.000      10.000  319.268    20.5     46
    0.005   False   10.000  27.000      10.000  151.535    21.9     47
    0.004   False   10.000  23.000      10.000   45.351      21     48


Thursday 2014 October 9
=======================

Summary
-------
* Code resolution checks and debugging (mainly `fglists.dat`)
* Submitted CRESST annual progress report snippet


Resolution notebook modifications
---------------------------------

Simplified some functions/framework for varying model computation (no need to
edit SNR parameters, just supply kwargs to fit functions).  I change the
bounding parameter values in notebook to make things work (using new fglists
has changed some FWHMs a good deal!).

To allow for larger tables, I change xex/fex in Fortran code to accommodate up
to 200 entries.  Also change disttab accordingly.

Full model code modifications
-----------------------------

To do computations, I modified the code to accept different filenames for
`fglists.dat`.  This requires us to re-read the table file on every calculation
for both Python code and f2py compiled Fortran.

    FullEfflength_port.py
        fefflen arg: fgfname (sent to FullEfflength_mod.f)
    FullEfflength_mod.f
        readfglists arg: fname
    models.py
        width_cont **kwarg: fgfname='fglists.dat' passed to fmp.fefflen

Next I pass through several internal integral resolution numbers...

    FullEfflength_mod.f
        distrmlt1, distrpohl, distrmgt1 args: itmax, inmax
        (distrmlt1 doesn't have inmax)
        distr args: itmax, inmax
        Fullefflengthsub args: itmax, inmax
    FullEfflength_port.py
        fefflen arg: itmax, inmax


Sources of model calculation error, in e- distribution table
------------------------------------------------------------

I first ran some tests using Sean's default `fglists.dat` with SN 1006
parameter bounds (not addressing damping, yet).

The effects to consider:
* Sampling (both number of samples, and sample spacing)
  It seems sensible to sample more finely near peak (most important
  contribution and more curvature in spectrum).  But, interpolation in the
  numerical calculation takes place in e- distribution, emissivity, and
  elsewhere.  So it doesn't seem obvious that, e.g., extremes of e- energies
  (low/high) are necessarily unimportant.

* Errors in Pacholczyk values vs scipy computed values

  Comparing hand-input Pacholczyk table values, with scipy computed values, the
  errors in Pacholczyk values are generally 1% or so (probably consistent with
  round-off error), but between x=3 and x=5 errors are as much as ~4-5%.

  For reference, `scipy.special.kv` is drawn from the AMOS Fortran library, in
  `zbesk.f` (last revised 1989?..) so I don't know if there are better routines
  available.  But it ought to do.

* Extending table values (to larger range of xex vlaues)

  Appending entries at new y values (larger or smaller) changes FWHMs by a LOT.
  I.e., 10% or more.  This is a huge problem.  Explore this by adding entries
  at both ends of Sean's tabulated values.


Error analysis for default fglists.dat
--------------------------------------
Scribbled notes.

`scipy.special.kv` table gives max error 0.22%, typical errors 0.05% (median),
0.07% (mean).

Log-space sampling (w/ `scipy` computed values) gives max error 0.9%, typical
errors 0.1-0.2% (median), 0.2% (mean).  So error due to resampling is important
(at 35 table entries).  This holds even if we compare resampled FWHMs to
hand-sampled FWHMs computed w/ scipy values (i.e., scipy errors are
comparatively small).

If we extend the number of entries, at the high frequency tail (1e1 out to
5e1):
* add `1.2e1`: max error is 10% (!), median 0.03-0.06%, mean 1.3%.  Maximum
  change is ~10% in all bands
* add `2e1` (alone, no value at 1.2e1): must change rminarc to 100 to fit one
  set of FWHMs.  The error on that one is huge.  Max error is 42%, median error
  is 0.1-0.3%, mean error is 5-6%.  It's just one point that blows up so badly
  and skews a lot of the errors actually, for mu=0 and eta2 = 1000 (both B0
  values).  The next largest FWHM errors here are about 5%, for mu=0 and eta2 =
  100.
* add all entries in Pacholczyk between 1e1 and 5e1: max error is now 73%
  (again, only at that one set of points).  Otherwise the max error is now
  about 1-5%, again at the mu=0/eta2=100 set of points.  The median error is
  0.08-0.2%, mean error 5% must be heavily skewed.
* Add entries 1e-5, 5e-5 (rest of table same up to 1e1): maximum error is
  0.0025%, Median error is ~1e-12, mean is ~1e-6.  Okay, awesome, there's no
  need to go below 1e-4.  Not surprising.

NOW, performing a similar analysis w/ doubled # of entries between 1e-4 and 1e1
(35 to 68/70)

* Doubling number of entries (35 to 68) with hand selection gives max error
  0.3%, typical errors 0.06% (median), 0.08% (mean).
* Doubling number of entries (35 to 68) using `scipy` computed values w/ manual
  sampling gives max error 0.2%, typical errors 0.03% (median), 0.05% (mean).
* Doubling number of entries (35 to 70) with log-sampling and `scipy` values
  gives max error 0.6%, typical error 0.05% (median), 0.1% (mean).  What if we
  compare to hand-sampled entries + `scipy`?  Max error becomes 0.5%, typical
  errors about the same or a bit smaller (0.05% median, 0.08% mean).

So again, log-sampling vs. hand-sampling incurs the largest fractional change
(but now the fractional change is a bit smaller than before).  It's not obvious
which one is more accurate, though.

__Conclusion:__ use a smart manual sampling, use more entries, and use scipy
values.  The big questions are 1. error due to adding/cutting pts at tail end
of distribution, and 2. error due to poor sampling of distribution.

New e- table sampling in `fglists_mod.dat`
------------------------------------------

Created new table with 89 entries, see the file for notes / code to regenerate.

Okay, I keep changing the sampling around as I run these resolution tests.
Key observations:
* xex values between 1e1 and 1e2 matter a LOT to the error, both in being
  present and in their sampling density

e- integral resolutions
-----------------------
For mu < 1, setting itmax below around 907 (default 1000) causes the
integration to fail miserably -- distr returns nothing but nans.

Investigating bug now...


Friday 2014 October 10
======================

Summary
-------
* Code resolution error checks cont. (internal integrals etc.)
* Apply new best settings.  Confirm Tycho, Tycho damping, SN 1006 errors
  generally 1% at max.


New resolution numbers / etc
----------------------------
After lots of resolution testing, and timing checks:

    irmax = 200     (from 400)
    itmax = 200     (from 1000)
    inmax = 50      (from 100)
    irhomax = 2000  (from 10000)

along with new `fglists_mod.dat`.  Also changed Tycho e- spectrum index `s` to
2.61 (for radio spectral index 0.58).

Looking at error results for damping -- error from e- distr integrals does
increase in damped calculation!  Which would explain why Sean chose itmax =
5000 in the damped code.

e- integral resolutions, cont.
------------------------------

Friday morning, caught the damned bug.  Turns out that at the end of the
integral, the rescaled variable n is supposed to be zero but can be like -1e-14
or 1e-14 or something.  In the negative case, this kills the integrand which
becomes NaN, thereby killing everything else!
I checked edge cases for the rest of the integrals.  I believe this shouldn't
happen elsewhere (other integrals have different limits etc.).

Now passing through `irhomax` (from Python port code) through as a parameter as
well, in `FullEfflength_port.fefflen` and `models.width_cont`.

### Timing check

Now we are using itmax=200, inmax=50, irhomax=2000, vs. 1000/100/10000 before.
We are using `fglists_mod.dat` with 2.5x more e- energy entries.
Running kernprof on the FWHM calculation (not considering adaptive calculation,
which would just be 2x longer), the time breakdown is now:

Tycho, 5 energy bands: 1.02 seconds

    77% time to fullmodel.distr
    7% time to emistab = emisx(...)
    4% time to emisgrid/intensity interpolation/integration
    11% time to FWHM computation

FWHM computation breakdown:
    15% time to spopt.minimize_scalar (find intensity max)
    42% time to spopt.bisect (right threshold)
    41% time to spopt.bisect (left threshold)

So, all looks good.


B-damping tables (set-up)
-------------------------

Set up generating script (see relevant script for details), ready to go.
Fixing shock speed vs to default.  Call by hand for a range of ab values.
The simple model initial guesses are pretty bad, but the FWHMs are relatively
insensitive to large values of B0 -- so that is good!

Modified table generating code to not get stuck when we're getting unable to
find FWHMs (and code returns 240 arcsec / 900 arcsec / remnant radius).
B0-stepping code detects when a FWHM is at the SNR's shock radius, indicating
error, and only fills in values between valid FWHMs.

Quick time estimate:
* 1 mu, 151 eta2, 50 B0 (ask for 30 points)
* 12 `a_b` values
* 2 sec/call (advective soln takes 1 sec/call)
* Run in "parallel", so divide by 3

Total time is: `1*151*50*12*2/3` = 60400 seconds or 17 hours.
This should be a slight overestimate.

`a_b` values I use are:

    0.5,    0.05,   0.04,   0.03,
    0.02,   0.01,   0.009,  0.008,
    0.007,  0.006,  0.005,  0.004

But, feed damping numbers into scripts as:

    python tab_tycho_damping_141010.py  0.5   0.01   0.005
    python tab_tycho_damping_141010.py  0.05  0.009  0.004
    python tab_tycho_damping_141010.py  0.04  0.008  0.007
    python tab_tycho_damping_141010.py  0.03  0.02   0.006

to try to balance the load of calculations (I have gone diagonally in this
table).  Started around 8:14pm, expect to finish by 1pm tomorrow.
So I can come in the afternoon and set it up to actually run fits.


Tycho downstream spectra
------------------------

I forgot to link downstream spectra to backgrounds, woops.  Fixed, and reran
fits to absorbed power law.

Approaches to fitting lines: 1. excise line energies, 2. fit lines with
Gaussians.

Modified `spec_fit.py` to fit/excise Sulfur line as well.
* Energies excised are 1.7-2 keV, 2.3-2.6 keV.
* Fixed Gaussian line fitting to properly keep previous freeze/thaw cycle
  parameters
* Changed parameter outputting (to npz/json/whatever) to avoid collision
  between models w/ two same components (keying of dict changed)
Reran ALL Tycho regions-5 spectra

For Tycho regions-5, I don't compute eqwidth for region 19.  The downstream
region is so tiny that it's useless for fitting.  Wait until next set of
regions (When I manually excise FWHMs) and then this region will become okay.

Added new table with fits/excisions (more formatting/patch-ups tbd, but good
enough for now).


Sunday 2014 October 12
======================

* Set up comprehensive magnetic damping fits
* Table output code for magnetic damping fits + mE values from best fits

Pretabulated damped FWHM tables
-------------------------------

General remarks: finished around 8-9am Saturday morning.  Lots of problems have
arisen.  Small glitch with deleting the tee'd standard output/error, but that
shouldn't matter so I haven't looked into it.  Some logs are weird (have
multiple tables' logs), but all the output is saved at least.

Observation: as damping gets stronger (ab ~ 0.03 to 0.01, for Tycho), we start
to see smaller B0 values being accepted.   Not only that, but at large eta2
values we see an increasing range of B0 values giving reasonable FWHMs.  Both
larger and smaller than the range of B0's at smaller eta2.


### Issue: calculations without damping enabled

When the initial guess for B0 is bad and the code computes a new value and
tries again, I forgot to pass the kwargs along (bug now fixed).  So many
calculations had no damping!  I reviewed the terminal output and found:

    ab = 0.006, eta2 = 83.18 to 100 are invalid (3 values)
    ab = 0.005, eta2 = 33.11 to 100 are invalid (13 values)
    ab = 0.004, eta2 = 7.24 to 100 are invalid (33 values)

For all larger values of ab, the initial guess for B0 was always accepted and
computations took place w/ damping.

I have modified these dictionaries BY HAND to remove the bad values, so that
free full-model fits don't attempt to use those values...

### Issue: insufficient sampling of small B0 / large FWHMs

B0 stepping overshoots badly, at large eta2 and moderate-strong damping.  The
FWHMs are insensitive to B0 changes at large B0, but at small B0 it starts to
matter much more.  But our linear stepping doesn't know this.  The following
eta2 ranges are strongly undersampled, having not enough small B0 values.

    ab = 0.01, eta2 = 75.86 to 100
    ab = 0.009, eta2 = 52.48 to 100
    ab = 0.008, eta2 = 33.11 to 100
    ab = 0.007, eta2 = 22.91 to 100
    ab = 0.006, eta2 = 13.18 to 100
    ab = 0.005, eta2 = 6.61 to 100
    ab = 0.004, eta2 = 0.79 to 100

This is a "continuous" problem; it affects computations at smaller ab / eta2 as
well.  In general, as eta2 increases the sampling of B0 gets coarser due to
this problem.  So a lot of tables don't cover parameter space that well.


Magnetic damping fits
---------------------

12 minutes for full/free fits at mu=1, for 22 regions... with a number of fits
running to massive eta2 (sigh).

With 12 ab values, and 6 eta2 values (1 free, 5 fixed) -- that's 12x12x6 = 14
hours... so it should finish before 9am.  Hopefully earlier because the fixed
eta2 fits should be faster.

We'll give it a shot, and see how it goes...
(status after approx 1 hr: 9/72 completed.  hopefully just 8 more hours, so
it'll be done in the morning for sure)

Set up code for printing out tables of magnetic damping fit results.
Also added code for to compute mE values for best loss-limited fits.


(Week 20) Monday 2014 October 13
================================

Summary
-------
* Assemble B-damping results


Magnetic damping results/tables
-------------------------------

Magnetic damping fits were quite fast!  Would have finished at around 11pm,
except the very last one (ab = 0.004) appears to have bugged up.  Generated
tables with damping results, `m_E` computed for best loss-limited and damped
fits, printed for tomorrow.  I went through results by hand and highlighted in
my tables.

Key observations:
* 13/20 regions are well-explained, w/ better chi-sqr values and a clear
  "valley" in chi-sqr space (eta2 vs. ab) for reasonable values of eta2
* 3/20 regions have damping fits w/ better chi-sqr values, but eta2 may be
  larger.  So we don't see the same trend
* 3/20 don't give damping fits w/ better chi-sqr values, but have damping fits
  within +1 of best chisqr.  Damping could be favored if eta2 is fixed to 1 (or
  similar).
* 1/20 seems to favor loss-limited model, by not showing any clear chi-sqr
  space trend (Region 18)

Basically, the damping model is able to explain our observations just as well
as a loss-limited model.  There aren't any cases where the chi-sqr difference
is so strong that we can favor one model over another.

SN 1006 damping tables
----------------------

151 eta2 values, 50 B0 values (x1.2 gives maybe 70), 12 ab values.
I estimate 1 second/call (3 bands vs 5) = 30.2 hr, so hopefully should be done
by tomorrow.


Tuesday 2014 October 14
=======================

Mostly assessing results, reading/thinking/math-ing today.  Not super
productive...

Meeting on damping results
--------------------------

* Spectra (lines / excision fits) -- just go ahead and state it in words.
* Results -- show the whole slew of results for a few regions, the best and
  worst.  Show that no matter what, the models will fit fine.  Plot data with
  errors and best fit curves.  At the end, give some kind of parameter range of
  best fits.

We have, in a sense, a more useful (certainly more correct) result.  Not as
useful or as satisfying a conclusion as we had thought, but still useful.

SN 1006 damping fits
--------------------
FWHM tables finished by this morning.  Set up same set of 72x damped model fits
for SN1006 (started around 8:45am).

Sent results around to Rob/Brian, after meeting.

Similar to Tycho -- Filaments 1, 2, 5 can be well-fit w/ a damped model.
Filaments 3, 4 (sub-Bohm in Sean's paper) aren't well fit w/ damping.  They
require loss-limited model + negligible diffusion to match the extreme drop-off
in rim widths.

But, again, we don't know if the individual region fits might differ...

Reviewing e- distr solutions (Lerche & Schlickeiser 1980)
---------------------------------------------------------

Walked through a little bit of Pacholczyk on plasma processes (not really
directly helpful), and went over the Green's function derivation by Lerche &
Schlickeiser.  This is basically reproduced verbatim by Rettig & Pohl.

The spatial dependence of `D(z) B^2(z)` mucks things up when we try to take a
Fourier transform to simplify.  The assumption is that

    \alpha = 1 + \frac{B_0 - B_{\mt{min}}}{B_\mt{min}} e^{-z/a_b}

is constant, to obtain analytic solutions.  This does become roughly constant
a few scale lengths behind the shock, when the B-field spatial dependence is
very weak.  But we are most interested in e- transport at and immediately
behind the shock, and I see no sign that we can neglect the exponential
variation there.

How can we circumvent this assumption, obtain an approximate solution, or show
that it doesn't matter?  Perturbation analysis, finite difference/element
solution, ???  I don't know the best way to go about this.

Wednesday 2014 October 15
=========================

Write everything up.  Material/questions for Brian, Sean, Steve soon.

srcutlog local model for XSPEC from Brian. Doesn't seem to initialize properly
at the moment, come back to this later...


Thursday 2014 October 16
========================

Summary
-------
* Some work on manuscript
* Attempt to fit with damping length `a_b` free

Paper updates
-------------
Integrate comments/markup from October 3 draft (last iteration that got sent
around to Rob/Brian).

Integrate comments on spectra fits (line fits / excision) into text.

Magnetic damping free fits
--------------------------
Modify model code to allow fitting on ab.
When calling `Fitter.get_fit`, add another kwarg `ab_fit` and be sure to enable
damping in `models.width_cont`.  Also modified fit saving code to handle `ab`.
Manual error bound code is NOT tested with ab.  Use at your own risk...

Tested out -- code runs into error saving.  But, it looks fairly useless
anyways.  The code is unable to explore chi-sqr space and stays very close to
the initial ab value (0.005, in this test).

Downstream spectra fits
-----------------------
Checking a few of the lower quality chi-squared values

* Region 5, redchi with 2 lines = 1.60
  Basically just a slew of lines.  I'm able to fit:
  0.84 keV (O/Ne mush?), 1.35 keV (Mg XI),
  1.86 keV (Si XIII, He alpha), 2.02 keV (Si Ly alpha), 2.22 keV (Si He beta),
  2.44 keV (S XV), 2.85 keV (S Ly alpha/He beta mush?),
  3.1 keV (Ar He alpha), 6.45 keV (Fe XXV whatever)
  Strongest line, other than Si/S He alpha, is that 3.1 keV Argon line.
  but generally just a lot going on
* Region 20 -- same problem, just a whole rigmarole of lines that have to be
  fit.  This one is even harder to fit than Region 5.


Random thoughts on B-damping (also, radio rims)
-----------------------------------------------

* Simple (catastrophic dump) equation can be formulated with magnetic damping
  too, entering via the `+ f/\tsynch` term (cooling timescale changes,
  lengthens as we move downstream).

  Naive attempts to solve ODE with `scipy.integrate.odeint` did not work;
  solutions blow up.  My _guess_ is that the solution is being dominated by an
  exponential-blowup-like solution (akin to exponential solutions for 2nd order
  ODEs, think modified Bessel functions or w/e).  My naive boundary conditions
  (f = 1, f' = 0 at shock) are most likely wrong, at least the f'=0 part.

  Better approach: try solving as a 2-point BVP, using shooting/relaxation
  methods.  We don't know what the interior boundary value is.  e- distr should
  decrease, but no justification for it to be zero at center (and there could
  be an age limit on propagation too -- ignoring all ejecta interaction etc).
  Or, w/ shooting method, maybe the increasing exponential will be forced to
  near zero and there will be effectively only one free parameter,
  normalization set at the shock.  So we wouldn't have to worry about the
  interior.  But try it first...

* Damping solution to transport equation could give better agreement with
  radio rims?  From punching in numbers briefly -- I think the transport model
  setup can reproduce these little mini-rims on top of a smush of radiation.

  Cassam-Chenai obtained radio profiles with rims, but the absolute intensities
  did not match data (trying to match radio / X-ray data concurrently).  The
  fall-off in radio peaks is also a bit steeper than some I've seen, so our
  models could be useful?

  Idea would be to simultaneously fit X-ray and radio profiles -- instead of
  just fitting FWHMs, fit all data together.  Should be more robust too.

  Reynoso et al. (1997) used VLA data at 1.4 GHz from 1994, 1995 (~1 arcsec
  resolution).  John Hewitt's VLA proposal stated 5 GHz in A configuration
  (max basline -> sub-arcsec resolution!)

  1 GHz     4.136 x 10^-9 keV
  1.5 GHz   6.2 x 10^-9 keV
  5 GHz     2.07 x 10^-8 keV

  Since e- energy losses are small, small differences in energy shouldn't
  matter much; and, the diffusion coefficient will be pretty small.
  Generate plots for range of parameters and see how they look.


Friday 2014 October 17
======================

Summary
-------
* Plot modeled intensity profiles (X-ray and radio)
* Reinstall heasoft (srcutlog works, PyXSpec no longer needs 32-bit Python!)

### Short drop-in to Brian's office

* parameter numbers/bounds -- yeah, eyeballing should be okay.  Look at the
numbers Monday and see.

* `srcutlog`: try sending spectra to Brian and he can compare fit #s real fast.
  Addendum: I reinstalled heasoft from source and it works `-_-`

Plots for paper restructure
---------------------------

Add region numbers to plots of spectra, profiles (to help keep straight)

* Region 1 -- no 0.7-1 keV profile, B-damping just barely better, chisqr ~ 25
* Region 10 -- B-damping barely worse, chisqr ~ 1 (redchi ~ 0.4)
* Region 14 -- only region where B-damping is a lot better than loss-limited
               chisqr ~ 10-20
* Region 16 -- B-damping just barely better, chisqr ~ 4
* Region 18 -- B-damping disfavored, chisqr ~ 10
* Region 19 -- B-damping just barely worse, chisqr ~ 90

How about -- walk through pipeline with Regions 1, 14, 18.

Modeled intensity profiles
--------------------------

Added flag to `FullEfflength_port.fefflen`, `models.width_cont`.
Breaks fitting functionality, of course!

Pass kwarg to `Fitter.width_full` or `models.full_width` and it will spit out
intensity data, r-grid values alongside FWHM output.  Need to also set
`rminarc`, `irad_adapt`, `iradmax` appropriately to get good profiles.

* Nice plots of tycho data with model fits (subplot if you need)
  (done)
* Nice plots of SN 1006 data with model fits
* RUN fits to individual SN 1006 regions, see what happens
* Run fits for SN 1006 AND Tycho with e- cutoff disabled, see what happens.
* Add text for Sean/Steve...
* srcutlog work.

* ODE solution to B-damping from simple model
  (tackled)

* MAKE profiles of radio/X-ray profiles for various damping parameters
  1. from our best fit numbers,
  2. from some "standard"/"nice" values
  3. to produce the funny radio shapes
  Done -- main observation is that strong B field amplification AND damping
  produces very sharp radio rims; weaker amplification and damping gives rise
  to shapes closer to what we might expect


Saturday 2014 October 18
========================

Summary
-------
* SN 1006 plots, fits w/ no e- cutoff
* Attempts to solve simple model ODE with B-damping

SN 1006 stuff
-------------

SN 1006 plots -- generated plots w/ loss-limited + damped predictions, profiles

Run fits without e- cut-off energy, to see if that accounts for the
expectation of constant rims (`m_E ~ -0.1`).
Result: box errors in many places for small damping lengths (`a_b`).
Essentially our initial guesses give FWHMs far too large, since without a
cut-off energy we get tons more electrons.  Results may not be reliable... but
investigate this later.

Naive IVP ODE solution to catastrophic dump equation
----------------------------------------------------
This one, with f'(0) = 0, gives reasonable numbers and goes straight to zero at
around `x = 40*ab`.  But, not sure if the shape is right.  It is clearly mixing
both the decaying and blowing-up solution, and I suppose there's no way to
discard the latter in our numeric solution except by fortuitous choice of
boundary/initial values...

Possible approaches considered:
* change of variables (results in only one non-constant coefficient, but it's
  specified in terms of a Lambert W / product log function)
* Series solution (not well explored).  Same issue of needing good choice of
  IV/BVs, using the Runge-Kutta or whatever is probably simpler here
* Laplace transform -- results in ugly convolution-like integral in
  non-constant coefficients.

But I haven't explored these that deeply.  Anyways here's the code snippet for
scipy's ODE solver.  See also `code/models/models.py`, since I suspect there's
a good chance no one will ever go through these notes in the future.

    from __future__ import division

    import numpy as np
    import scipy as sp
    from scipy.integrate import odeint
    import matplotlib.pyplot as plt

    CD = 2.083e19   # c/(3e)
    b = 1.57e-3     # Synchrotron cooling constant

    eta = 1  # Regular eta, not eta2
    mu = 1

    v = 1.25e8  # shock veloc, 1/4 of 5e8 = 5000 km/s
    E = 36.5    # e- energy for 1 keV photon at B = 100 \muG
    Bmin = 5e-6
    B0 = 30e-6
    ab = 1e17   # Damping length ~1% of Tycho's shock radius

    rs_tycho = 240 * np.pi/180. /3600. * 3.0 * 1e3 * 3.085678e16 * 1e2
    print rs_tycho

    def B(x):
        """Damped B-field"""
        return Bmin + (B0 - Bmin) * np.exp(-x/ab)

    def Bprime(x):
        """Spatial derivative of damped B-field"""
        return -1/ab * (B0 - Bmin) * np.exp(-x/ab)

    def D(x):
        """Diffusion coefficient"""
        return eta * CD * E**mu / B(x)

    def odefunc(vec, x):
        """Function for ODE derivatives at x """
        f, g = vec
        c1 = - b* np.power(B(x),2) * E / D(x)
        c2 = v/D(x) + Bprime(x) / B(x)
        return g, c1*f + c2*g

    ic = [1, 0]  # f(x=0) = 1, f'(x=0) = 0

    xvec = np.linspace(0,100*ab,1000)
    svec = odeint(odefunc, ic, xvec); fvec = svec[:,0]
    plt.plot(xvec,fvec); plt.ylim(0,1); plt.show()


(Week 21) Monday 2014 October 20
================================

Summary
-------
* Damping analysis plots/tables for Tycho/SN1006; paper updates
* Email to Steve/Sean on damping results
* Install CASA for EVLA data (symlinks in `/usr/bin/`)
* New Tycho tables with (1) more `a_b` values, and (2) Bmin = 1e-6

Morning meeting
---------------
* Email to Steve/Sean today (done), need their feedback (make sure we're doing
  it right, close the hole by having Sean confirm this)
* Discussion! Write it up (`-_-` should have been done last week but it's okay)
* ODE solving -- probably not going to yield much in that alley.
  Check literature, in case solution is already obtained (I think not, though)
* VLA data -- Jack gone this week, but can look through manuals and bombard him
  with questions next Monday!  Drag him into this project... (poor guy)

* Brian: for data plots, show the best region, w/ consistent decreasing trend,
  that we can get.
* Any consistent trends within filaments -- comparable B fields or eta2 values?
  Try fixing eta2 = 1 and get best fits for moderate damping, or with either of
  loss-limited/damped.

* Rob is out Thursday/Friday-ish.  Next Monday Brian out in the morning, meet
  after lunch.

Stuff to Sean & Steve (Rob/Brian)
-------------------------------
(Brian sent email ~2pm today)

Attached are plots of width-energy data and model predictions, for the
following two cases (SN 1006 figure has legend).  Damping length is fixed (ab =
length / shock radius) and injected e- spectra have exponential cutoffs.
Solid black curve is best damped fit with mu = 1, dotted/dashed curves are
loss-limited with varying mu.

SN 1006, Filament 1
    Best damped fit: ab = 0.006, eta2 = 0.04, B0 =  34.6 muG, chisqr = 0.15;
    best loss-limited fit: eta2 = 2.59, B0 = 102.5 muG, chisqr = 0.14
    (ab = 0.006 is 5.4" in SN 1006)

Tycho, Region 16.
    Best damped fit: ab = 0.004, eta2 = 0.22, B0 = 341 muG, chisqr = 3.58;
    best loss-limited fit: eta2 = 3.82, B0 = 583 muG, chisqr = 4.25.
    (ab = 0.004 is 0.96" in Tycho)

Big questions: 1. can Sean replicate our results (how was `m_E ~ -0.1`
determined), and 2. are the transport equations correct?  Specific points:
* Does an exponentially varying diffusion coefficient make sense?
  A diffusion gradient term `-dD/dx * df/dx`, from `-d/dx (D(x) df/dx))`
  becomes important, and, where damping is strong, increasing diffusion
  coefficient can even cause rims to become _narrower_!
* How can we assess the validity of the assumption that `D(x) B^2(x)` is
  constant? (needed to derive equations (14-16) in Sean's paper)
* Could other assumptions be undermined by the presence of damping (plane
  shock without curvature, constant downstream velocity, ???)
* If equations are right -- can we safely apply the transport model to radio?

Note to Sean: I am using slightly different resolutions + a larger
`fglists.dat` (~2x resolution for single e- emissivity); SN 1006 fits use 3
data points.

SN 1006, compare results with/without e- cutoff
-----------------------------------------------

Summary: yes, e- cutoff makes it harder to get good fits (although should be
noted that we are now comparing loss-limited and damped fits, both without e-
cutoff).  Only Filament 5 was able to find a better damped fit.

Values of `m_E` to compare between loss-limited fit, damped fits with/without
cutoff energy.  All evaluated at mu = 1.  Damped fits are best fits for ab less
than 0.01 (where we tried free fits at 0.01, 0.009, 0.008, ... 0.004).

            Loss-limited fit    Damped fit, cutoff        Damped, no cutoff
            ----------------  -----------------------  -----------------------
    Region  mE(1-2) mE(2-3)   ab      mE(1-2) mE(2-3)  ab      mE(1-2) mE(2-3)
    ---------------------------------------------------------------------------
    1       -0.35   -0.31     0.006   -0.34   -0.30    0.005   -0.27   -0.26
    2       -0.49   -0.48     0.004   -0.69   -0.42    0.010   -0.24   -0.27
    3       -0.51   -0.51     0.010   -0.63   -0.41    0.006   -0.59   -0.38
    4       -0.51   -0.50     0.008   -0.65   -0.47    0.004   -0.61   -0.34
    5       -0.19   -0.18     0.010   -0.20   -0.18    0.006   -0.26   -0.17

Funny how steep the `m_E` values are for the best damped fits, even in the
case without damping!  Curious.

In general, the damped fits with no e- cutoff have eta2 < 1 (sub-Bohm).  A lot
of the free fits failed pretty badly, but I think the identified "best fits"
are reasonably representative of the actual best fits (looking at fits with
eta2 fixed, they're not any better).


Sean's reply and more plots/work
--------------------------------

Generate plots showing variation in profiles at low energy (overlay of profiles
in loss-limited and damped cases; width-energy relation to 0.1 keV).

What plots do we have so far:

    energywidth         data + model predictions, width vs. energy
    energywidth-lowE    0.1 to 5 keV, mu=1 only
    prfs                profiles for all energy bands fitted
    prfs-radio          radio and X-ray bands
    prfs-modelcomp      kevs = np.array([0.05, 0.1, 0.2, 0.5, 1., 2., 4.])

These are somewhat time consuming to run (especially the energywidth plots).
Regenerated all Tycho, SN1006 plots w/ new annotations etc as needed.

Generate new table with Bmin = 1e-6 (other parameters same)

Paper plots/tables layout
-------------------------

Following some meeting discussion, re-generate all relevant tables and plots.

* FWHM table
* Best loss-limited fits (mu varying)
* Best magnetic damping fits (ab varying) (state DOFs in footnotes)
* Plots of fitted profiles (use dashed lines for loss-limited fit)
* Table of best damped + loss-limited fits for all regions, side by side
* Table of best damped + loss-limited fits as above, with eta2 = 1 fixed!

Now throwing out:
* table of loss-limited full model fits for all regions, all mu values
* table of alternate spectral fits
* table of filament-averaged parameter values
Some of these tables are useful to draw numbers from, and can/should be
regenerated on the fly.  But they won't make it into the final paper.

Planning to throw out simple model fits completely, at this point.  But just
comment out in text for now, don't discard.

Updated tables/plots to use Tycho Region 16.  Added data to compare
loss-limited / damped fits side-by-side (incl. with eta2 = 1)


Tuesday 2014 October 21
=======================

Summary
-------
* Review plots, address Sean's suggestions
* Paper writing
* CASA tutorial (EVLA observation of 3C391) meddling

Working with damping results
----------------------------
Ran numbers/stuff by Brian.  Seems like we're at a bit of a question mark
point.  Edited plots of different model FWHMs to show half-max line,
illustrating where the measurement is being taken.  Replied to Sean's email
with plots (fits with Bmin = 2e-6 pending).

Hold off on generating Regions 6 w/ culled regions.  If we end up fitting
profiles, rather than FWHMs, that will change our strategy entirely.

Doesn't seem obvious to me there's any trend within filament numbers.

### Ramblings

Some observations: the best fits with small magnetic field are the ones that
give us strong plateaus behind the rim, + give "spurious" energy dependence due
to measurement being taken near the base of the rim.

How about fits where damped rims are thin and don't show a plateau, so the
energy dependence is not so suspect?
E.g., region 14 w/ the best decrease in chi-squared value w/ a damped model.
Best fit damped length is 0.005 or 1.2 arcsec, rim widths are 2-3 arcsec.
The rims show practically the same energy dependence at our X-ray energies, as
the loss-limited model, with eta2 = 100 and eta2 = 3 (not unreasonable; if we
set eta2 = 1 then B0 = 50 is needed with still pretty good chi-sqr).

Answer: the width-energy dependence is pretty weak, so no wonder -- it does
just fine.

Where energy dependence is weak, no one cares.  Magnetic damping works.
Where (best-fit-curve) energy dependence is stronger, magnetic damping goes
weird / tricksome.  BUT this is the only way to get smaller B fields.
We note too that the magnetic damping in "good" cases still gives
pretty strong magnetic fields, if smaller than in the loss-limited case.
(this is, roughly, the train of thought Parizot followed)

That is interesting -- that in these cases with more mild energy dependence, we
still see sizable fields in the damped model.
... think more about this tomorrow.

### PDE assumptions...

Poked around with full model PDE work, seeing if anything could simplify.  Can
rewrite in terms of Lambert W function, somewhat more simply (just passing
through the integral transform in x/z-coordinate), but still looks intractable.

I really don't know what to say about the `D(x) * B^2(x)` constant assumption.

Thought: I am assuming the cut-off energy is set by `B = B_0`.  Should it be
set at `B = B_min` instead?  Since behind the shock the diffusion coefficient
is highest at low `B`... don't know if that makes sense.  And Bohm limit is
supposed to be a minimum, we really need more constraints on this somehow.

### Workarounds

Spectral fits downstream?  To try to capture how we see higher plateaus behind
damped rims...



Damped fits with smaller Bmin
-----------------------------

Tables done, set up fits to start running overnight


Wednesday 2014 October 22
=========================

Summary
-------
* Parameter space exploration (bad behavior for moderate energy dependence)
* Paper -- update results/text, organize/start writing discussion (now
  accepting that B-damping is okay)


Parameter space exploration
---------------------------

Damping results with Bmin = 2e-6.  Not much difference is seen.

Generated grid of plots to help visualize parameter space, in X-ray and radio
both.  Only care about X-ray for now but I'm keen to see what happens in the
radio (shapes much more diverse, may be more helpful in constraining parameter
space...).

Regenerated plots for paper too (following style/layout used for this grid)

Drop-in with Brian
------------------

### Damping stuff

Ran plots of parameter space and thoughts by Brian.  Idea: just X out all the
plots that look non-sensical.  Just explain what we see; it would be useful if
we could say something about what fits are favored, what damping fits are okay
or are not okay. 

Paper first, please -- it will help everyone see which holes need to be
filled in, and assess how to do so.  Remember that we are focusing on an
observational result -- leave more detailed/correct modeling to others!

Get something as complete as possible -- clean up the little notes and other
details.

### Radio emission in SNRs

Q: any other sources for radio emission?  E.g., in the remnant interior.<br>
A: you still need GeV electrons for radio synchrotron emission.  The only thing
that might compete with our rims is acceleration at the reverse shock.


New Tycho tables
----------------

Generating new tabs with the new default resolutions for both
loss-limited and damped cases.  Then we'll be ready to go with new numbers
tomorrow...


Thursday 2014 October 23
========================

Paper!

Error on damping fits
---------------------

Trial error computation on the best fits on Regions 1, 16, with `a_b` fixed
Changed adaptive computation algorithm -- now step by % change in chisqr, up
until within +/- 1 of threshold.  Prevent it from getting stuck far from
threshold.

Also now accepts brentq kwargs, so user can set relative/absolute tolerance on
errors (same for both eta2/B0, should be ok given that our errors are quite
large -- tolerance in bounds should be like 1% of (bound - true\_val)).

Errors (saved as test output data somewhere as well)
(bounds now good to 0.01%)

* Region 1, ab = 0.008

    Best fit: B0 = 24.55, eta2 = 2.01, mu = 1
    Error bounds:
        B0 = [19.123, 24.55, 31.751 muG]
        eta2 = [0.48, 2.01, 11.04]

    Best fit: B0 = 19.21, eta2 = 11.18, mu = 2
    Error bounds:
        B0 = [17.45, 19.21, 21.48]
        eta2 = [4.33, 11.18, 30.63]

* Region 16, ab = 0.004

    Best fit: B0 = 340.54, eta2 = 0.220, mu = 1
    Error bounds:
        B0 = [326.4, 340.54, 355.1]
        eta2 = [0.112, 0.220, 0.381]

    Best fit: B0 = 318.6, eta2 = 0.235, mu = 2
    Error bounds:
        B0 = [300.5, 318.6, 336.1]
        eta2 = [0.107, 0.235, 0.461]

A big problem is that we're not varying `a_b`.  But, the best fits at other
values of `a_b` may give a rough estimate of the spread in `B_0`, and that
spread is somewhat large.


Friday 2014 October 24
======================

Summary
-------
* Paper draft with discussion/conclusions to Rob and Brian
* Plot of best damped fit `m_E` values (takes a little time to generate)
* CASA tutorial (for VLA data)


(Week 22) Monday 2014 October 27
================================

Summary
-------
* Review / markup paper.  Minor fixes to manuscript
* New Tycho FWHM tables finished
* Inspect Tycho best fits w/ smaller damping lengths (ab = 0.003, 0.002)
* CASA tutorial on VLA 3C391 data, up to CLEANing step

Paper markup
============

Small fixes.  I totally forgot paragraph on SN 1006 results (and how they would
change).

New Tycho tables, regions-6, tweaks to plots/whatever
=====================================================

### `m_E` dichotomy due to poor FWHM measurement

Generated better version of plot showing `m_E` values for best damped fits, as
function of energy.  Unequivocally shows that damped fits behave like Sean
explained, but only the spurious / edge cases blow up like crazy!

Regions 1, 2, 4, 5, 12, 18 (Tycho regions-5) all exhibit this blow up.

IF I repeat all of this, now with new fits at ab=0.003, 0.002, will other
regions do this as well?  Most likely, yes.  So let's try it.

But, remember, it's a spurious effect tied to us just picking that best fit, it
doesn't mean that other good (non-spurious) fits don't exist!
Is there a systematic way to sieve out these bad fits?
(e.g., require that damped rim is present at 0.1 keV or 0.01 keV?)
(but that also imposes a limitation, on downstream emission essentially)
(kind of arbitrary selection criteria)

### New Tycho tables

New tables (started Wednesday October 22) have Bmin = 5e-6, more data pts
(asked for 50 B0 values), more ab values (0.002, 0.003).  Damping tables all
completed on 2014 October 24.  Loss limited tables all completed on October 23.

Checked error logs, all clear.

Now using tables to generate new fits for ab = 0.003, 0.002 for regions-5, with
Bmin=5e-6.  Region 18 seems to be getting stuck on the free eta2/B0 fit.
I am wondering if smaller damping values could give smaller B-field values for
the Tycho regions w/ weaker energy dependence.


Tuesday 2014 October 28
=======================

Summary
-------
* Meddling with CASA CLEAN (done)
* Meeting w/Rob, Brian on paper

Running CLEAN, put together some plots/tables for SN 1006 to look at this more
closely...

Meeting notes
-------------

* About ready to send around to Steve/Sean, minor fixes... how to explain SN
  1006 results (discussion of dichotomy in fits looks okay)
* CRESST retreat is a go (good dry run for AAS)
* Leave for Dec 5 is fine
* No meeting next week, Nov. 3 (both Brian, Rob on plane at that time)
* Nov. 10 Rob most likely out; Nov. 17, 24 should be around (24 is thanksgiving
  week already!)

Next steps?
* Kepler -- finish passing through the pipeline
* VLA stuff ... maybe in brian/jack's paper(s) somewhere
  Remains low priority but especially pending feedback from Steve/Sean
* Poster next week!


Paper updates - plots and fixes
-------------------------------
Generated 4x5 subplot of energy-width curves (lot of twiddling plot settings,
tickmarks, etc)

Refactored plotting code, not optimal but much cleaner now.

Incorporate comments from Brian and my own notes


Wednesday 2014 October 29
=========================

Summary
-------
* Paper sent to all collaborators
* Start pipeline for Tycho regions-6, loss-limited fits

Misc: download/install astropy, APLpy to make nice Tycho image (publication
quality)

Paper
-----

Finished incorporating all comments, add text, plots, figures on SN 1006
Now describing the dichotomy in energy-dependent fits by "weak-field" damping
and "strong-field" damping.


Tycho Regions 6
---------------

Almost identical to regions-5, but I cut out regions 21, 22 and apply cuts to
FWHMs (blacklist the known bad quality ones -- though I'm afraid we could be
throwing out good data, it is consistent with our procedure and we mention this
smapling bias in the discussion).

Generated:
* `regions-6.reg`
* `regions-6.physreg`
* `regions-6-box.reg`  (for display purposes)
* `regions-6-az.txt`

Profiles generated with prefixes `prfs`, `prf-cts`, `prf-proc`
Processed profiles with fit domain cuts:
* `profiles/prf-proc_fit_cuts.dat, profiles/prf-proc_fit_cuts.npz`

Applied FWHM blacklist with `code/profiles/profile-fit-tycho.ipynb`
* `fwhms-subbkg/fwhms_*` (and `fwhms_spec_cuts.[dat,npz]`)
* `fwhms-subbkg/plots/prfs_*.pdf`

Run loss-limited and damped model fits with new tables:
* `tables/Tycho_gen_2014-10-23_grid_6-151-20_vs-*.pkl`
* `tables/Tycho_gen_2014-10-24_grid_1-151-50_ab-*.pkl`
Shock velocities interpolated from:
* `data-tycho/tycho_velocs.txt` (Williams et al. 2013), d = 3 kpc assumed
One-particle synchrotron emissivity taken from
* `code/models/fglists_mod.dat`
All default settings used (copy-paste from `models.py` and `snr_catalog.py`):

    def width_cont(params, kevs, snr, verbose=True, rminarc=None, icut=None,
        irmax=None, iradmax=None, ixmax=None, irad_adapt=True, irad_adapt_f=1.2,
        idamp=False, damp_ab=0.05, damp_bmin=5.0e-6, fgfname='fglists_mod.dat',
        itmax=200, inmax=50, irhomax=2000, get_prfs=False):

    tycho.dkpc = 3.0
    tycho.rsarc = 240
    tycho.s = 2.16 # e- spectrum index, 2.16 = 2*0.58 + 1
                   # 0.58 = radio spectral index from Sun et al. (2011)

    tycho.vs = 3.6e8 * tycho.dkpc/2.3  # Overriden by interpolated vs values
    tycho.cratio = 4.0

    tycho.icut = True
    tycho.rminarc = 20
    tycho.irmax = 200
    tycho.iradmax = 100
    tycho.ixmax = 500

    tycho.par_init = {'mu': 1.0, 'eta2': 1.0, 'B0': 300e-6}
    tycho.par_lims = {'mu': (0., 2.),
                      'eta2': (1e-16, 1e5),  # model code div by zero on eta2=0
                      'B0': (1e-6, 1e-2)}

No kwargs passed to `scipy.leastsq`

Run with the commands (after usual initialization):

    mu_vals = [0., 1./3, 0.5, 1, 1.5, 2.]
    outroot = TYCHO_OUTDIR + 'full-man_err'

    f_data = mdsp.build_dataf('full', conf_intv=0.683)

    dview.push(dict(f_data=f_data, mu_vals=mu_vals, outroot=outroot))
    plists = dview.map_async(lambda args: mdsp.generate_fits(f_data,
                args[0], args[1], mu_vals, outroot, save=True), gen_tycho())

Started around 8:52 pm today.


Thursday 2014 October 30
========================

Summary
-------
* srcutlog calculations of eta2 (equation for `nu_cut` fixed)
* synchrotron constant check
* generate regions-6 spectra

regions-6 spectra set up to generate overnight (up regions first)

srcutlog
--------
Fitting functionality added to `spectra/spec_fit.py` (fit type 3)

Add cell in plotter-prfs-spec to tabulate srcutlog fit parameters
(quickly get sense of fit quality and break frequency, compare to plain fits)

New scripts:
* code to interpolate shock velocities from azimuth angles, Brian's table
  (note: must supply distance scaling by hand, not hard-coded)
* code to parse break frequencies from `spec_fit.py` output, and compute
  table of eta2 values

Error in existing `nu_cut` formula -- off by a factor of 100 microGauss.
Also combined two constant terms w/ same exponent.  Correct eqn reads:

    \nu_{\mt{cut}} = c_m
        \left( 13.3 \unit{erg} \right)^{4/(1+\mu)}
        \left( 100 \muG \right)
        \left( 2657 \unit{erg^2} \right)^{-(1-\mu)/(1+\mu)}
        \left( \frac{v_s}{10^8 \unit{cm/s}} \right)^{4/(1+\mu)}
        \left( \eta_2 \right) ^{-2/(1+\mu)}

Now computed values are reasonable, and match those of Parizot (2006)
The only thing left to do is to run fits with eta2 values fixed, but regions-6
is still being computed.

Synchrotron derivations
-----------------------

Looking through synchrotron constants again and got caught up in details.
Finally figured out that Sean's `c_m \approx 0.29 c_1` (see Condon and Ransom).
When averaging single electron power over solid angle, the factor of
`\sin\vartheta` (where cursive theta is pitch angle) just makes life suck.

But happily all the constants check out nicely.  I realized that
`\langle \sin^2 \vartheta \rangle = 2/3` and so all stated constant values
check out.  The only question is of the synchrotron emissivity, and the
delta function invocation when converting electron cut-off energy to
synchrotron break frequency.

Little factors don't matter much here.  Back to important things...

APLpy figure generation
-----------------------

Installed `pyregion` from `pip` so APLpy can parse DS9 regions.
`pyregion` can't parse projections, must be boxes.

Notes:
* maybe disable TeX usage in matplotlibrc!  Too slow otherwise
* run `ipython --pylab` or change the backend to OS X or whatever, to get the
  figure to display properly

Seems like we must generate RGB images (png, FITS in DS9 first,
then pass it on to APLpy

    import aplpy
    f = aplpy.FITSFigure('rgbimg.fits')  # Saved from DS9
    f.show_rgb('rgbtest.png')            # Saved from DS9
    f.show_regions('regions-6-box.reg')
    def ra2deg(h,m,s):
        return (h+m/60.+s/3600.)/24.*360.
    f.recenter(ra2deg(0,25,18), 64.14, width=10./60, height=10./60)


Friday 2014 October 31
======================

Summary
-------
* Regions-6 fits and spectra
* APLpy image for paper
* Kepler full model resolution check, spectrum fits, FWHM table

Misc things from Brian
----------------------

Flux question -- stupid question, it was right in the `merge_obs`
documentation...  units are photons/cm^2/sec.
Convenient unit conversion
[form](http://www.stsci.edu/hst/nicmos/tools/conversion_form.html)
Okay to experiment around a bit, waiting on our collaborators anyways.

Brian and Rob are in Japan... send poster by end of day Tuesday (their Weds
morning), they can look it and send back in the same day.

Tycho regions-6 (cont)
----------------------

Regions-6 rim and downstream spectra merged, linked to backgrounds.
Backgrounds taken from `data-tycho/bkg-2`
Rim spectrum fits: phabs x powerlaw, phabs x srcutlog; alpha = 0.58
Downstream fits: phabx x powerlaw plain, lines excised, lines fit.
All look good (enough counts in all sections, rim or downstream)

Loss-limited fits with full errors are done.
Saved tables of interpolated shock velocities, srcutlog-derived eta2 estimates
Set up fits with eta2 values from srcutlog fixed (very quick, only took ~20
minutes).

APLpy figure (cont.)
--------------------

Pulled out some files/notes on figure generation that were created/kept in
`data-tycho/regions-4`, to use for CRESST/AAS posters.  Easier to find now.
The images/settings should be independent of the regions I'm using.

Kepler: 5kpc, shock radius 1.9 arcmin, shock veloc 4700 km/s assuming 5kpc at
ear, radio index -0.64.  Brian suggested just averaging the measurements
between Katsuda/Vink (Katsuda: 4700 km/s from Regions 4,5; Vink: 5250 km/s from
some kind of sector-wise full fitting routine, though he treated the SE ear
filament separately).  Average is 4975 km/s which I just call 5000 km/s for
simplicity.

Kepler resolutions, first FWHM tables
=====================================

Note: last Kepler work was October 7 -- I generated regions-3 and sent it
through the pipeline, simply needed to run full model fits.

BUT, I don't have spectrum fits yet!  The spectra are linked to background
files (`data-kepler/bkg-1`).  I run fits for rim and downstream spectra (same
procedure as for Tycho regions-6) with alpha = 0.64.

Current setup: regions-3 has 7 region selections on Kepler's SE ear, and one
sitting on the more mangled NW ear.  Satoru's proper motion paper didn't sample
the mangled region, and Vink's paper does it by sector -- so that area's motion
gets averaged out with the slower NW stuff behind it.

So for simplicity I will assume 5000 km/s throughout, no need for tables at
multiple velocities.

Resolution -- see notebook.  All looks good for Kepler



(Week 23) Monday 2014 November 3
================================

