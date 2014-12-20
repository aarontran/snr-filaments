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
* Week 22 - (10/27) Tidy, augment, send manuscript. FWHM fits w/ srcutlog eta2
            Run Tycho regions-6, Kepler regions-3 through data pipeline
* Week 23 - (11/3) Emails (param space), spectral variation, VLA tutorials,
            CRESST retreat and poster
* Week 24 - (11/10) Inspect radio data, more spectral variation, VLA RFI
            tutorial

* Week 25 - (11/17) ..., ..., first fits to full xray/radio profiles
* Week 26 - (11/24) Radio profile shape mapping
            (Thanksgiving, out Weds-Fri [work remotely Weds])

* Week 27 - (12/1)
            (out Friday)
* Week 28 - (12/8) Finish eyeball fits, send paper around, ..., ..., poster
                   minor loose ends
* Week 29 - (12/15) Send poster around, ...
* Week 30 - (12/22, 12/23 only)

* Week 31 - (12/29)
            (at home)
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
-------------------------------------

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

Summary
-------
* Data pipelines (Kepler loss-lim, damped fits)
* CASA tutorial, small Tycho EVLA observation (TDEM0020, Pannuti + students)
* Plots of FWHM vs. E in parameter space (Steve's email)

Sean and Steve's comments
-------------------------

Started assembling plots (to address Steve's questions/comments)

Kepler table
------------

Kepler tables finished early Sunday morning.  Log observations:
ab = 0.003, 0.002 seem to have had some trouble w/ getting B0 points.
ab = 0.004 got 15 pts at eta2 = 100, ab = 0.005 got 41 pts.  So those are okay

Set up for Kepler loss-limited and damping fits, running now.
At very small ab values the fits don't appear to move very much, and I'm not
sure if that's due to the fit or due to the bad starting parameters (from FWHM
table).


Radio data meddling
-------------------

Continuing CASA 3c391 tutorial with multiscale cleaning
Python crashed unexpectedly partway through iterations... (was running pretty
slowly to begin with)

> unhandled exception: Unrecognized map name '' 
> libc++abi.dylib: terminate called throwing an exception
> Caught an ABORT  signal. Please wait ...Abort trap: 6

Downloading very recent TDEM0020 data from September 2014
(`TDEM0020.sb29665410.eb29703711.56924.434018796295.ms`)
Demo Tycho observation by T. Pannuti and Morehead State students


Tuesday 2014 November 4
=======================

Summary
-------
* Draft CRESST poster
* Damped fits with mu = [0, 1./3, ..., 2], Tycho regions-6
* Varied plots of parameter space for Steve/Sean email
* Kepler regions-3 damping fits, inspect some plots/tables

Most time today generating/twiddling plots, writing email response... need to
work more efficiently.

Sean/Steve's comments (cont)
----------------------------

New set of parameter grids (profiles), plots of FWHM-energy dependence w/
varying parameters, plots of `a_b / FWHM` (no particular trend/constant
behavior, we observe a good range).

Sent email with plots, discussion of damped fits with variable mu

CRESST poster draft
-------------------

Modified and updated NASA poster appropriately.  Presenting equivocal
results, but AAS poster should have more.

Tycho regions-6, Kepler regions-3 pipeline
------------------------------------------

Tycho damped fits with eta2 = 1 fixed, multiple mu values (Steve) (done)
Running slew of damped fits for regions-6 overnight
Then run loss-limited fits with eta2 = 1 fixed.

Kepler regions-3 damped fits complete.  Starting to inspect plots and tables...
results look comparable to Tycho.

    Region & $\eta_2$ (-) & $B_0$ ($\mu$G) & $\chi^2$
           & $\eta_2$ (-) & $B_0$ ($\mu$G) & $\chi^2$ & $a_b$ \\
    \midrule
    1 & 752.07 & 1373.1 &  1.95 &     0.08 & 204.9 &  1.62 & 0.006 \\
    2 &   2.70 &  446.2 &  4.80 &     2.70 &  22.9 &  4.45 & 0.003 \\
    3 & 599.67 & 1251.9 &  8.63 &     0.01 &  46.4 &  6.06 & 0.003 \\
    4 &   2.40 &  403.5 &  3.73 &     0.00 & 290.3 &  3.30 & 0.008 \\
    5 &   0.16 &  337.0 & 32.41 &  2021.65 &  10.1 & 30.18 & 0.006 \\
    6 & 531.97 & 1351.3 & 60.92 &     0.00 & 247.8 & 58.78 & 0.005 \\
    7 & 667.60 &  912.7 &  3.61 &     0.00 &  87.3 &  0.24 & 0.007 \\
    8 &   0.00 &  266.0 &  1.73 & 15941.01 &   7.6 &  0.02 & 0.008 \\

Regions 1, 3, 4, 6 (looks like crap), 7 are acceptably fit with a "strong"
damping model.  Regions 2, 5, 8 give "hybrid" fits to reproduce strong energy
dependence (borne out by plot of `m_E` vs. energy, hitlist is 2, 5, 8).
Funnily all the damped fits are better, though.


Wednesday 2014 November 5
=========================

Summary
-------
* Inspect rim/downstream photon index variation
* XSPEC fakeit to test effect of curved spectrum
* Review physics of thermal (ionized) line emission (see new md file)
* New spectrum fits for 2-7 keV only -- larger change in photon indices

Tycho regions-6 damped fits completed this morning

Spectral variation between downstream/rim
-----------------------------------------

### First look at spectrum fits

Following Sean's suggestion, try examining differences between rim and
plateau energy spectra.  Rettig & Pohl (2012) give predictions of spectral
index for loss-limited + strong damping cases (varies smoothly over soft
X-rays).

Set-up: update spectrum table code (remove hacky workaround to ignore Region
19, for Tycho regions-5). Updating a few tables to Regions 6 in manuscript
(inconsistently) so I can look at spectrum fit parameters quickly.

The average spectral indices are:
* thin rim: 2.8043
* downstream, excised: 2.8406
* downstream, line fit: 2.8427

The differences between indices are 0.036 (excise) and 0.038 (line fit).
Not bad!  Compared to Rettig & Pohl, this could be consistent with either
loss-limited rims or weak damping.  Hard to say, though.

### Check results for power-law fit to smoothly steepening thing

Because Rettig/Pohl predict the spectral index to smoothly increase (spectrum
gets steeper), I wonder how our fits are affected by this.

I attempt to fakeit in XSPEC with a 2-break power law.  I download the ACIS-I
responses from [CXC](http://cxc.harvard.edu/caldb/prop_plan/imaging/), using
the Cycle 16 on-aimpoint files (`acisi_aimpt_cy16.rmf, acisi_aimpt_cy16.arf`).
I require a normalization of 1e-4 (similar to our spectra) and 750 ks exposure.
Here is a loss-limited model:

__Loss-limited fakeit__

    ========================================================================
    Model phabs<1>*bkn2pow<2> Source No.: 1   Active/Off
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.600000     +/-  0.0          
       2    2   bkn2pow    PhoIndx1            2.25000      +/-  0.0          
       3    2   bkn2pow    BreakE1    keV      1.00000      +/-  0.0          
       4    2   bkn2pow    PhoIndx2            2.40000      +/-  0.0          
       5    2   bkn2pow    BreakE2    keV      3.00000      +/-  0.0          
       6    2   bkn2pow    PhoIndx3            2.60000      +/-  0.0          
       7    2   bkn2pow    norm                1.00000E-04  +/-  0.0          
    ________________________________________________________________________

If I fit this to an absorbed power law, best fits are:

    ========================================================================
    Model phabs<1>*powerlaw<2> Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.758752     +/-  2.33789E-02  
       2    2   powerlaw   PhoIndex            2.73490      +/-  3.85094E-02  
       3    2   powerlaw   norm                1.35289E-04  +/-  5.91311E-06  
    ________________________________________________________________________

    ========================================================================
    Model phabs<1>*powerlaw<2> Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.600000     frozen
       2    2   powerlaw   PhoIndex            2.50952      +/-  1.95295E-02  
       3    2   powerlaw   norm                1.01600E-04  +/-  1.60192E-06  
    ________________________________________________________________________

__Damped fakeit__

Okay, looks plausible.  What about if I consider the damped filament?
The spectral indices are smaller overall.

    ========================================================================
    Model phabs<1>*bkn2pow<2> Source No.: 1   Active/Off
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.600000     +/-  0.0          
       2    2   bkn2pow    PhoIndx1            2.10000      +/-  0.0          
       3    2   bkn2pow    BreakE1    keV      1.00000      +/-  0.0          
       4    2   bkn2pow    PhoIndx2            2.25000      +/-  0.0          
       5    2   bkn2pow    BreakE2    keV      3.00000      +/-  0.0          
       6    2   bkn2pow    PhoIndx3            2.45000      +/-  0.0          
       7    2   bkn2pow    norm                1.00000E-04  +/-  0.0          
    ________________________________________________________________________

Fit faked data to absorbed powerlaw now.

    ========================================================================
    Model phabs<1>*powerlaw<2> Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.746972     +/-  2.23844E-02  
       2    2   powerlaw   PhoIndex            2.60887      +/-  3.58063E-02  
       3    2   powerlaw   norm                1.40123E-04  +/-  5.79448E-06  
    ________________________________________________________________________

    ========================================================================
    Model phabs<1>*powerlaw<2> Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.600000     frozen
       2    2   powerlaw   PhoIndex            2.40638      +/-  1.81683E-02  
       3    2   powerlaw   norm                1.08007E-04  +/-  1.62775E-06  
    ________________________________________________________________________

__fakeit conclusions__

The overall, very rough shift on filament spectrum appears discernable as a
shift in best fit PhoIndex of ~ 0.1.
 
XSPEC fitting code (spectral variation)
---------------------------------------

### Procedure

Fit spectra between 2.0 and 7.0 keV.  Freeze nH=0.7 (see Hayato et al. (2010),
Cassam-Chenai (2007?)).  First try this automatically (with excised or line
fits).  Then, try manual fits to account for extra lines if needed (in rim and
downstream spectra).

In Region 5 of Tycho regions-6, besides Si/S He alpha the next strongest line
is 3.1 keV Ar (from notes, 2014 October 16, toying around in XSPEC).
Just to be aware of.

### Work/results

Refactored PyXSPEC code somewhat, to make fit information dumping/etc more
practical.  Confirmed that fits w/ cleaned code look identical.

Added spin-off script (`code/spectra/spec_var_fit.py`) to generate fits over
domain 2-7 keV.  Then quickly parse results in IPython.
The results appear even stronger now!  Average photon indices are:

* thin filament:    2.92 (std=0.15)
* downstr excise:   3.23 (std=0.14)
* downstr fit:      3.17 (std=0.19)

Average differences are thus 0.31 (excise) and 0.26 (fit).
Stdev of differences are 0.17 (excise), 0.22 (fit).

So this appears consistent with weak to moderate damping, compared to the
loss-limited model.  It's certainly not rigorous, but I think it's compelling.

### Manual fitting (Ar He alpha and more)

Inspect automatic line fit plots -- which regions could benefit most from
manual fitting?

* Region 2 -- some bump at ~ 3 keV (possible S He beta, Ar He alpha)
* Region 5 -- clear 3.1 keV Ar line, maybe 3.9 keV Ca
* Region 6 -- maybe 3.9 keV Ca
* Region 7 -- maybe 3.1 keV Ar
* Region 8 -- 2.9 keV S He beta? 3.1 keV Ar? maybe Ca..
* Region 9 -- maybe 3.9 Ca
* Region 10 -- try 2.9 keV S He beta
* Region 11 -- try S, Ar
* Region 13 -- Si He beta (2.2 keV), Ly beta (2.38 keV), Ar (3.1 keV)
* Region 14 -- Ar (3.1 keV, 3.7 keV?)
* Region 15 -- try 3.9 Ca
* Region 19 -- Si Ly beta (2.38 keV)
* Region 20 -- S He beta, Ar He alpha

Summary -- most are worth a shot.  Strong contenders differ, but S He beta, Ar
He alpha, maybe Ar He beta / Ca He alpha.

Try Region 5 first.  The automatic fits gave:

    Region 5
        rim, rim-excise, rim-fit: 3.047, 3.342, 3.318
        3.047, 0.2954, 0.2710

which is actually pretty average (probably large # of counts help).

Attempting to fit an Si Ly alpha line makes a big difference!...
Current fit has:
- 2.21 keV Si He beta (bad attempt)
- 2.44 keV S He alpha
- 3.12 keV Ar He alpha
PhoIndex changes from 3.275 to 3.125 with/without 2.00 keV Si Ly alpha.
Huge difference!  Because the line is right at the edge of the fit, it
sort of tips the falloff around too easily...

But, I'm not sure if the fit is even physically reasonable.
Here we have:

    1    1   phabs      nH         10^22    0.700000     frozen
    2    2   powerlaw   PhoIndex            3.12517      +/-  8.62544E-02  
    3    2   powerlaw   norm                2.27262E-04  +/-  2.55635E-05  
    4    3   gaussian   LineE      keV      2.00000      frozen
    5    3   gaussian   Sigma      keV      4.30906E-02  +/-  5.00851E-02  
    6    3   gaussian   norm                1.07679E-06  +/-  6.04577E-07  

vs. the fit without:
 
    1    1   phabs      nH         10^22    0.700000     frozen
    2    2   powerlaw   PhoIndex            3.27552      +/-  6.34065E-02  
    3    2   powerlaw   norm                2.79041E-04  +/-  2.13870E-05   

The fitted gaussian has with 0.04 keV and has normalization 1/3rd that of the
Sulfur line.  But, we expect Ly alpha lines to be narrower (not sure how much
broadening to expect.  Looking at Hayato et al., the line is very weak and
might well be covered over by the Si He alpha line's spread.  BUT, we are
looking right behind the shock -- could the near-shock Si be more strongly
ionized?

    ========================================================================
    Model Model Component  Parameter  Unit     Value
     par  comp
       1    1   phabs      nH         10^22    0.700000     frozen
       2    2   powerlaw   PhoIndex            3.32061      +/-  5.86091E-02  
       3    2   powerlaw   norm                2.97017E-04  +/-  2.04703E-05  
       4    3   gaussian   LineE      keV      2.44056      +/-  7.06256E-03  
       5    3   gaussian   Sigma      keV      3.28679E-02  +/-  1.23684E-02  
       6    3   gaussian   norm                2.43370E-06  +/-  2.60479E-07  
       7    4   gaussian   LineE      keV      3.11902      +/-  0.179961     
       8    4   gaussian   Sigma      keV      1.00787E-03  +/-  7.46441E-02  
       9    4   gaussian   norm                5.18801E-07  +/-  1.07458E-07  
    ________________________________________________________________________

Verdict -- trying to fit Si lines near edge is fishy (can cause changes of ~0.1
in photon index), but adding stuff in the middle is probably not so bad.

I might well trust the excised fits more at this point, which actually gives a
stronger spectrum difference.  Let me try one more, now excising the Argon line
too to be very very safe (cut 3.0 to 3.1).

Excising S line (2.3-2.6 keV) only: mean PhoIndex = 3.2311
Excising S, Ar lines (2.3-2.6, 3.0-3.2 keV): mean PhoIndex = 3.2316

Summary: negligible...?


2014 Thursday November 6
========================

Summary
-------
* Ran few more spectral variation fits (testing..), inspect data
* Plan for full model spectra to go w/ best fits
* Print CRESST poster
* Finished 3C391 VLA tutorial


Spectral variation summary
--------------------------

Here are a few examples.  These are photon indices for phabs x powerlaw fits,
for, in order:
* rim spectrum (2-7 keV)
* downstream spectrum (2-7 keV) with
  - 2.3-2.6 keV excised
  - 2.3-2.6, 3.0-3.2 keV excised
  - 2.45 keV line fit
* downstream spectrum (2.6-7 keV)

And the differences between the downstream and rim spectra

    Region 1
        rim; excise, excise_ar, line, tail: [ 2.845  3.219  3.221  3.172  3.17 ]
        diffs: [ 0.374  0.376  0.327  0.325]
    Region 5
        rim; excise, excise_ar, line, tail: [ 3.047  3.342  3.345  3.318  3.255]
        diffs: [ 0.295  0.299  0.271  0.208]
    Region 16
        rim; excise, excise_ar, line, tail: [ 2.92   3.359  3.36   3.318  3.067]
        diffs: [ 0.439  0.44   0.398  0.147]

For comparison, if I fit the rim spectrum between 2.6-7 keV ("tail"), I get
photon indices of 2.889, 3.124, 2.875 for Regions 1, 5, 16.

Now looking at just the averages (fit over 2-7 keV unless otherwise stated):

    Fit type           Photon index        Diff fr rim      Diff fr rim (2.6-7)
    ---------------------------------------------------------------------------
    rim                2.918 (std 0.151)   -                -
    rim, 2.6-7 keV     2.936 (std 0.152)   -                -
    down, excise S     3.231 (std 0.137)   0.31 (std 0.17)  0.30 (std 0.17)
    down, excise S/Ar  3.233 (std 0.138)   0.32 (std 0.17)  0.30 (std 0.17)
    down, fit S line   3.173 (std 0.185)   0.26 (std 0.22)  0.24 (std 0.22)
    down, 2.6-7 keV    3.166 (std 0.180)   0.25 (std 0.19)  0.23 (std 0.18)

These are averaged photon indices, averaged differences (same as difference of
averages, but we can get a std-dev from indiv. differences).

Differences are just downstream index - rim index.  So we do see that the rim
spectra are harder in most cases -- and to be precise, I haven't even gone
through and munged the fits by hand or anything

Here are the full data for reference 

               Rim indices           Downstream indices*
             ----------------   ------------------------------
    Region  2-7keV  2.6-7keV   cut S   cut S/Ar  fit S   tail
    ------  -------------------------------------------------
    1       2.84    2.89       3.22    3.22      3.17    3.17
    2       2.98    3.05       3.44    3.45      3.41    3.32
    3       2.83    2.92       3.25    3.24      3.25    3.22
    4       2.95    3.04       3.2     3.19      3.18    3.25
    5       3.05    3.12       3.34    3.35      3.32    3.25
    6       2.93    3.14       3.12    3.12      3.12    2.97
    7       2.95    2.99       3.22    3.22      3.05    3.47
    8       2.79    2.76       2.86    2.86      2.86    2.94
    9       2.92    2.95       3.13    3.14      2.59    3.04
    10      2.76    2.92       3.02    3.02      3.01    3.09
    11      2.82    2.69       3.09    3.1       3.09    2.66
    12      2.51    2.5        3.22    3.22      3.21    3.03
    13      2.85    2.87       3.35    3.35      3.31    3.34
    14      3.08    2.97       3.35    3.36      3.24    3.34
    15      3.03    2.9        3.26    3.27      3.26    3.3 
    16      2.92    2.87       3.36    3.36      3.32    3.07
    17      2.97    2.97       3.33    3.34      3.34    3.19
    18      3.3     3.13       3.16    3.16      3.14    3.07
    19      2.86    2.91       3.27    3.27      3.25    3.27
    20      3.03    3.1        3.41    3.41      3.32    3.32

    *Fit on domain 2-7 keV, except tail is 2.6-7 keV

Inspecting the data by hand -- only in Region 18 is the downstream spectrum
consistently __harder__ than the rim spectrum (spectral index smaller by 0.14
to 0.23), everywhere else we see clear softening (one or two isolated places
where the downstream spectral index might be smaller than upstream, depending
on fit type being considered).  The effect is not as dramatic if we compare to
the 2.6-7 keV rim fit.

This is not predicted in either the loss-limited or damped cases, so something
may be up with Region 18, but it's not obvious to me what exactly is up, based
on profiles/spectra/region location.

How can we best present this tabular data??  Right now it's probably easiest to
just state in text and compare to the model predictions/output.  Might want to
show rim/downstream spectra overlaid on each other.

Modeling spectral variation
---------------------------

How to use our model to discuss/investigate this?
As Sean attempted before -- we can extract intensity based on best fits.

The downstream rim extraction region is very variable, depended on fits.
In model, go from the 0.7 keV FWHM out back by one more FWHM distance.
(so the regions should have the same length.

So, take our fits -- extract the radial intensity profiles with fixed rminarc
at multiple, finely spaced energies.  This gives us keV/cm^2/s/keV...



Friday 2014 November 7
======================

Summary
-------
* CRESST retreat
* New (1994) VLA image from Steve (from D. Moffett?)
* Mucking around with test EVLA Tycho observation (TDEM0020, Sep 2014)
* Generate rest of Kepler tables for Monday -- very quick
* Generate fits (and all derived products) to assess Kepler spectral variation

Notes
-----

Comparing VLA and Chandra images, the expansion between 1994 and 2009 is very
visible.  Very cool.

Spent some time exploring Tycho TDEM0020 dataset (strong RFI everywhere),
pulled antenna corrections / basic data flagging / inspect RFI contamination.
All set to run through tutorial on automatic RFI excision.

Remark: Kepler data, when rerunning fits, seem to sometimes not find FWHM even
if repeating same calculation.  Adaptive FWHM blows up right at edge I guess...
Not precisely sure why this occurs (and, why this didn't occur in Tycho fits?)

For fun, run the spectral variation tests on Kepler very quickly.  Already have
it set up.  So we get a slew of "spectral variation" fits:

    python ~/snr-research/code/spectra/spec_var_fit.py regions3_up 0 plot_var/plt_var_up fit_var/fit_var_up -v
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_up 4 plot_var/plt_var_up_tail fit_var/fit_var_up_tail
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_down 0 plot_var/plt_var_down fit_var/fit_var_down
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_down 1 plot_var/plt_var_down_excise fit_var/fit_var_down_excise
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_down 2 plot_var/plt_var_down_excise_ar fit_var/fit_var_down_fit
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_down 3 plot_var/plt_var_down_fit fit_var/fit_var_down_fit
    python ~/snr-research/code/spectra/spec_var_fit.py regions3_down 4 plot_var/plt_var_down_tail fit_var/fit_var_down_tail

Average spectral indices are:
* 2.77 (rim, 2-7 keV)       std = 0.13
* 2.83 (rim, 2.6-7 keV)     std = 0.26
* 3.31 (down, excise S)     std = 0.25
* 3.33 (down, excise S/Ar)  std = 0.26
* 2.96 (down, fit line)     std = 0.55
* 3.11 (down, 2.6-7 keV)    std = 0.37

Average differences (for excise S, excise S/Ar, fit, 2.6-7 keV respectively):
(compared to rim spectrum fit over 2-7 keV)
* 0.54 (std 0.27) excise S
* 0.56 (std 0.28) excise S/Ar
* 0.19 (std 0.60) fit
* 0.34 (std 0.42) 2.6-7 keV

If we compare to the rim spectrum fit over 2.6-7 keV, subtract 0.06 from all
these numbers and ramp up the standard deviations slightly.


(Week 24) Monday 2014 November 10
=================================

Summary
-------
* Spectral variation computed from full model code

Misc: skimming VLA/CASA tutorials, threw Kepler tables into a (crappy) PDF
for easier review.  Didn't get quite enough done today...

Spectral variation from full model profiles
-------------------------------------------

(edit, 2014 November 12: the tabulated results are WRONG -- I did not properly
change the integration length to be over the 1 keV FWHM only!
Please see notes from November 12 for corrected tables)

Procedure:
* compute intensity profiles from best fit parameters over 25 log-spaced
  energies between 1.05-7 keV.
* Using the largest FWHM (at 1.05 keV), integrate intensities from shock
  to 1 FWHM behind shock (rim), and from 1 to 2 FWHM behind shock (plateau)
* Plot integrated rim/plateau intensities vs. energy to get energy spectra
* Fit spectra between 2-7 keV and back out spectral indices
* Subtract rim/plateau indices to get difference


Just wrangling data/presentation at the moment.
The spectral indices are very scattered.  But, interestingly, the spectral
differences are very small in the damping model.

Generally, I observe that the spectral index difference gets _larger_ at lower
energies for our best damped fits.  This strongly contrasts with the results
that Rettig/Pohl present.

Rettig/Pohl use w = 5 arcsec (assuming d = 2.5 kpc), which is very large -- the
plateau region is _well_ behind the rim.

### Spectral variation, best damped fits (`a_b < 0.01`), Tycho regions-6

Mean change is 0.035 (std = 0.030).  Maximum change is 0.088 (Region 12).
Perhaps unsurprisingly, the largest spectral index change is associated with
(1) weak-field damping, and (2) smallest damping length.

    Region	rim	    down	delta	ab	    B0	    eta2
    1	    -2.626	-2.672	0.046	0.008	25.4	1.63
    2	    -3.452	-3.540	0.088	0.003	16.8	19
    3	    -4.790	-4.840	0.051	0.003	14.1	184
    4	    -2.643	-2.691	0.047	0.004	25.2	1.72
    5	    -2.388	-2.459	0.071	0.003	27.6	0.487
    6	    -1.680	-1.688	0.0075	0.003	114.4	0.000356
    7	    -2.019	-2.021	0.0016	0.004	241.6	0.0172
    8	    -2.080	-2.097	0.016	0.003	52.2	0.152
    9	    -1.802	-1.841	0.039	0.002	46.8	0.0176
    10	    -2.271	-2.306	0.035	0.002	36.1	0.363
    11	    -2.216	-2.298	0.083	0.002	31.0	0.2
    12	    -2.020	-2.100	0.08	0.002	35.0	0.0518
    13	    -1.930	-1.927	-0.0032	0.004	243.4	0.00072
    14	    -1.790	-1.805	0.015	0.003	89.9	0.0119
    15	    -2.916	-2.935	0.019	0.003	31.3	4.48
    16	    -2.425	-2.416	-0.0096	0.004	340.5	0.22
    17	    -2.349	-2.350	0.001	0.003	98.4	0.463
    18	    -9.919	-9.977	0.058	0.006	8.7	4.1e+03
    19	    -1.847	-1.894	0.047	0.002	45.8	0.0241
    20	    -2.622	-2.632	0.0097	0.003	47.3	1.7

### Spectral variation, best damped fits w/eta2=1, Tycho regions-6

Mean change: 0.020 (std=0.025)
Max, min are 0.067, -0.008.  A lot more (barely) negative changes this time.

    Region  rim     down    delta       ab      B0      eta2
    1       -2.510  -2.550  0.04        0.008   27.7    1
    2       -2.534  -2.592  0.058       0.003   25.4    1
    3       -2.723  -2.714  -0.0083     0.010   406.5   1
    4       -2.513  -2.552  0.039       0.004   27.7    1
    5       -2.554  -2.621  0.067       0.003   24.5    1
    6       -2.459  -2.463  0.0033      0.004   68.0    1
    7       -2.597  -2.598  0.00079     0.006   290.9   1
    8       -2.489  -2.490  0.00069     0.005   137.8   1
    9       -2.461  -2.468  0.007       0.004   78.3    1
    10      -2.504  -2.544  0.04        0.002   29.1    1
    11      -2.689  -2.693  0.0033      0.010   316.8   1
    12      -2.541  -2.597  0.056       0.003   25.9    1
    13      -2.591  -2.591  -0.00054    0.006   257.7   1
    14      -2.508  -2.517  0.0093      0.004   59.8    1
    15      -2.514  -2.521  0.0064      0.003   47.8    1
    16      -2.722  -2.724  0.0019      0.006   434.1   1
    17      -2.578  -2.578  0.00062     0.004   232.1   1
    18      -2.602  -2.664  0.062       0.003   23.9    1
    19      -2.494  -2.494  -0.00023    0.004   77.6    1
    20      -2.494  -2.504  0.01        0.003   57.5    1


### Spectral variation, best loss-limited fits w/ eta2=1, Tycho regions-6

Mean change: 0.0005 (std = 0.006).  Max/min are 0.01, -0.01.
This all seems consistent with a net change of zero -- which sounds reasonable.

    Region	rim	    down	delta	    B0	    eta2
    1	    -2.776	-2.784	0.0079	    183.3	1
    2	    -2.778	-2.783	0.0058	    313.0	1
    3	    -2.778	-2.784	0.0055	    427.0	1
    4	    -2.780	-2.786	0.0058	    284.8	1
    5	    -2.783	-2.785	0.0027	    288.8	1
    6	    -2.784	-2.786	0.002	    410.9	1
    7	    -2.785	-2.783	-0.0016	    419.0	1
    8	    -2.785	-2.784	-0.00052	388.2	1
    9	    -2.786	-2.781	-0.0041	    415.1	1
    10	    -2.786	-2.782	-0.0043	    466.3	1
    11	    -2.789	-2.787	-0.0028	    355.6	1
    12	    -2.790	-2.787	-0.0036	    318.0	1
    13	    -2.792	-2.793	0.0006	    400.3	1
    14	    -2.811	-2.804	-0.0069	    383.9	1
    15	    -2.805	-2.815	0.01	    431.8	1
    16	    -2.807	-2.816	0.0088	    493.4	1
    17	    -2.810	-2.818	0.0085	    467.8	1
    18	    -2.812	-2.801	-0.012	    283.1	1
    19	    -2.803	-2.801	-0.0022	    401.4	1
    20	    -2.805	-2.794	-0.011	    463.6	1

But, anyways, we do see that the loss-limited spectra are even more invariant
than the damped ones.

### What have we found so far?

The measurements in Tycho suggested a typical change of 0.2 to 0.3 (spectral
index at rim ~ 2.9, downstream ~ 3.2).

The best model fit parameters are showing much smaller changes in spectral
index.  Damped fits see changes of around 0.02-0.03, loss-limited fits see
almost no change.

So something is going on here!  I don't know what, yet.

I erred earlier in "extracting" full model spectra from the FWHMs for each
energy -- e.g., at 4.5 keV I'd get spectra from rim to FWHM(4.5 keV), and
FWHM(4.5 keV) to 2x FWHM(4.5keV).  I changed this to use just FWHM(1.05 keV) as
the baseline.  But this made very little difference, except for giving more
sensible results in the loss-limited spectral variation.

I don't understand how Rettig/Pohl got their figures.

Tuesday 2014 November 11
========================

Summary
-------
* RFI removal from TDEM0020 (tutorial), running overnight (bad parameters...)
* Parameter range parsing for full model code, run overnight

Automated/manual RFI removal in CASA
------------------------------------

This is tricky and subjective.  See my notes in the TDEM0020 folder.

Mapping parameter space
-----------------------

Moved previous code for mapping parameter space, over to new notebook.
Now figuring out what to do...

What are the key qualitative behaviors that we know / have described?
* Smaller B0 -> widening rims (and, more likely to show hybrid behavior?)
* Smaller eta2 -> wider rims (small ab), or thinner rims (large ab)
* Smaller mu --> marginally thinner rims, but not important
* Smaller ab --> smaller rims

Now, the tricky part is getting past the zero-th/first-order effects.
* as eta2 -> 0, we get diminishing marginal returns.  This limits the max width
  of damped rims, or the minimum width / energy-dependence of loss-limited rims
* in strong damping, as B0 -> large, we see weakening effect.  The rims
  converge to a loss-limited width beneath the damping lengthscale at
  sufficiently high energy.  As B0 -> small, rim width blows up quickly.

Questions to answer:
* what combinations of parameters give rise to this "hybrid" model behavior?
  The key to differentiating is that they give comparable FWHMs at X-ray
  energies, but at lower photon energies (UV to visible) the FWHMs blow up.

  This would be a 2 step approach:
  1. find all parameters that give "reasonable" FWHMs at, e.g., 1 keV
  2. for all those parameters, compute FWHM at 0.01 keV.  If it blows up, it's
     not a "real" damping model.

* What range of ab values can reproduce a given profile, as Steve asked?


* Strong B-field, small ab value can set _upper_ bound on rim width
  (thick rims not possible; strong-field damping)
* Weak B-field, small ab value can set _upper_ bound on rim width
  (thick rims buried by plateau above 50%; "blowup/hybrid" weak-field damping)
* Weak B-field, moderate to large ab value cannot set _upper_ bound on rim
  width.  Rim will simply widen, but will have measurable FWHM
  (this encompasses the loss-limited case)

In general, for a given damping length, there won't be a _lower_ bound on rim
width -- just crank up magnetic field as far as possible.

Also -- it would be nice to distinguish between blowup (plateau above 50%) and
rims simply widening very far back.
 

Wednesday 2014 November 12
==========================

Summary
-------
* spectral variation, round 2 (debugging and validation)
* radio profiles, first look w/ plotting script

(Philae landing day)

Meeting with Brian
------------------

Brought material on Kepler, spectral variation, parameter space mapping, quick
ad hoc plots of radio and X-ray profiles together.

On modeling radio rims:
* Relic electrons in the radio -- see his RCW 86 paper, they can diffuse out
  and may persist for hundreds of years!
* Jack -- here Friday for sure...

* Spectral index variation -- yeah, interesting in that we explored, ran into
  trouble. get sean's feedback.  Show reproduced version of Rettig/pohl plots
  and explain that we can't use this?  Another paragraph for manuscript
* Radio spectral index variation -- too close together, probably can't get a
  good number out (large error bars).

* Radio -- draw attention to the fact that we don't see rims.  Ask for feedback
  explicitly on that.  Send plots around with rims/x-ray comparison
* Send SN 1006 stuff (just had telecon -- discussing results) (done)

Subsequently I compiled an email and sent radio profiles, better version of
parameter grid.  Now easier to pattern match by eye.


Spectral variation, round two (fixed, old tables bad)
-----------------------------------------------------

Caught a major bug -- I was not integrating properly over the same FWHM in each
region.  Fixed, changed our results somewhat.

Spent some time trying to reproduce Rettig/Pohl results.  No go on damping!
If I supply smaller magnetic fields (55 muG vs. 85) it looks much closer.  But,
the qualitative behavior is still wrong, and the 1 keV spectral index
difference is still wrong.  So something is clearly up.

### Best loss-limited fit, eta2=1 fixed

    Region	rim	    down	delta	ab	    B0	    eta2
    1	    -2.593	-2.638	0.045	nan	    183.3	1
    2	    -2.596	-2.641	0.045	nan	    313.0	1
    3	    -2.596	-2.642	0.045	nan	    427.0	1
    4	    -2.598	-2.643	0.045	nan	    284.8	1
    5	    -2.599	-2.644	0.045	nan	    288.8	1
    6	    -2.600	-2.645	0.045	nan	    410.9	1
    7	    -2.601	-2.646	0.045	nan	    419.0	1
    8	    -2.602	-2.647	0.045	nan	    388.2	1
    9	    -2.602	-2.648	0.045	nan	    415.1	1
    10	    -2.603	-2.648	0.045	nan	    466.3	1
    11	    -2.607	-2.652	0.045	nan	    355.6	1
    12	    -2.608	-2.653	0.045	nan	    318.0	1
    13	    -2.612	-2.657	0.045	nan	    400.3	1
    14	    -2.634	-2.678	0.044	nan	    383.9	1
    15	    -2.636	-2.680	0.044	nan	    431.8	1
    16	    -2.638	-2.683	0.044	nan	    493.4	1
    17	    -2.641	-2.685	0.044	nan	    467.8	1
    18	    -2.633	-2.677	0.044	nan	    283.1	1
    19	    -2.623	-2.668	0.045	nan	    401.4	1
    20	    -2.626	-2.671	0.045	nan	    463.6	1

    Mean 0.0449188299054
    Stdev 0.000446049746373
    Max 0.0454939453366
    Min 0.0441799959736

### Best damped fit, eta2=1 fixed

    Region	rim	    down	delta	ab	    B0	    eta2
    1	    -2.383	-2.450	0.067	0.008	27.7	1
    2	    -2.389	-2.466	0.077	0.003	25.4	1
    3	    -2.544	-2.594	0.05	0.010	406.5	1
    4	    -2.387	-2.454	0.068	0.004	27.7	1
    5	    -2.398	-2.481	0.082	0.003	24.5	1
    6	    -2.380	-2.414	0.034	0.004	68.0	1
    7	    -2.474	-2.513	0.039	0.006	290.9	1
    8	    -2.406	-2.437	0.031	0.005	137.8	1
    9	    -2.385	-2.418	0.033	0.004	78.3	1
    10	    -2.389	-2.453	0.063	0.002	29.1	1
    11	    -2.532	-2.579	0.047	0.010	316.8	1
    12	    -2.404	-2.477	0.074	0.003	25.9	1
    13	    -2.477	-2.515	0.037	0.006	257.7	1
    14	    -2.430	-2.465	0.035	0.004	59.8	1
    15	    -2.431	-2.471	0.04	0.003	47.8	1
    16	    -2.569	-2.615	0.046	0.006	434.1	1
    17	    -2.485	-2.518	0.033	0.004	232.1	1
    18	    -2.449	-2.530	0.081	0.003	23.9	1
    19	    -2.416	-2.449	0.033	0.004	77.6	1
    20	    -2.416	-2.452	0.036	0.003	57.5	1

    Mean 0.0504027841573
    Stdev 0.0179311849968
    Max 0.0823741798834
    Min 0.0314484727336


Thursday 2014 November 13
=========================

Summary
-------
* Radio tutorial on TDEM0200, trying to excise RFI more carefully
* Email chain discussion...

Some JVGR work

Email chain discussion
----------------------

Explanation for lack of radio rims in loss-limited model, from Sean:
we'd need to include self-similar solutions for B, v\_plasma, etc. behind the
shock.

Maybe: define metrics for various regimes
* When do we see a thin rim within 10 arcsec of shock, in the radio?
* When is the said thin rim just a thin bump (minimum behind shock about 80%)


Friday 2014 November 14
=======================

* Radio work (TDEM0200), excision and calibration etc
* JVGR work

Picked up Steve's conference proceeding paper (eds. Roger and Landecker)


Sunday 2014 November 16
=======================

Working on calibration for TDEM0020.  Not looking so great at the moment.
See my notes.  The calibration is very poor.


Monday 2014 November 17
=======================

Meeting deferred (Rob out, Brian Boston-bound early)
Primarily JVGR work


Tuesday 2014 November 18
========================

Summary
-------
* Clean up notes

Morning JVGR work

Short drop-in with Rob... not much new, but ran by some ideas/plots, and ideas
for mapping parameter space.  Monday, Rob may or may not be in, will know soon
(Rob in Cambridge, UK rest of week; Brian in Cambridge, MA).

Read/scan/save Steve's paper (Roger/Landecker) briefly.  It seems like Steve
explored the rim problem, accounting for some of these sphericity etc. effects
(and varying magnetic fields somewhat?).  But, you simply can't get rims as
thin as are seen in SN 1006.  Not directly helpful, but just a reference.

Fitting entire profiles
-----------------------

First iteration setup of code.  More tomorrow.

Considered and rejected iterative closest point type algorithms (iteratively
compute best rotations/translations).  Although I've done some of this before,
it's not worth it.

How to line up / align / transform profiles, intensities, etc?

  One idea -- use FWHMs to line up profiles, roughly (since we don't precisely
  know the shock location + precursor may be seen ahead of shock...)
  (consider the work of Warren 2005?)
  (morlino/caprioli 2012, cassam-chenai 2007 have done this to some extent)
  Another: match up profiles to peaks.  But where/how to cut-off front/back for
  fitting?...


Wednesday 2014 November 19
==========================

Summary
-------
* Finish profile fitting code
* Explore fit behavior / compare to data

Profile fitting
---------------

Code set-up and now working -- good number of bad quality fits.
Radio fitting is of course much faster because diffusion is unimportant, and
the advective solution is easy.

Current set up can fit multiple profiles w/ single scaling+translation.
But for now I just fit them separately, we have some issues to iron out first.
(different issues in profile fitting for each case)

Remarks/observations
--------------------

* Quality of best fit is generally very poor
* Very hard to depress plateau emission in full model, even w/ damping
* X-ray emission is always lower (after normalizing) than radio -- no surprise
* Modeled rise (esp. in X-ray) seems much sharper than measurements -
  precursor?
* Fits often not well constrained -- slowly iterating through a valley in
  chi-squared space by tuning B0, ab.

Need more robust indicators/measurements of shape!

Approaches and things...
* The best fit model cannot have emission ABOVE the plateau of our measurements
  (unphysical, unless our model is wrong)
* Can we fit behind peak to avoid precursor like stuff?
* Compute minimum emission behind shock
* Map out parameter space based on local minimum behind peak (and, presence or
  lack of a peak!)
* Maximum FWHM allowed for a given damping length.  I think we can get this.
* Self-similar Sedov-Taylor solutions to get slightly more correct downstream
  profiles (add a flag to trigger this, then compare with/without)
* Do we need to worry about smearing, deconvolution effects in radio image?
  I assume CLEAN already handles the beam PSF

EXAMPLE that fits NW rim (num=2) well in radio at bump, but falls off too fast
vs. observation (downstream/relic e- not considered in our model...)

        ab = 0.0165
        B0 = 381.447e-6
        amp = 3500
        r_trans = 0.87

        # X-ray
        amp=55
        r_trans = 1.5

Thursday 2014 November 20
=========================

JVGR

Friday 2014 November 21
=======================

JVGR

Saturday 2014 November 22
=========================

## Summary: pulling out profile shape parameters

Idea 1: constrain rims / damping based on downstream emission behind thin rim.
How far does the rim drop in radio?  From my sample regions (NNE, NW, WNW, SW),
I observe three cases for emission following initial rise at shock:

1. emission keeps rising, but more slowly (d^2I/dx^2 < 0) (loss-limited)
2. emission falls and keeps falling (weak damping)
   slow fall, to ~70% over 40 arcsec
3. emission falls, turns around (trough) and rises again (strong damping)
   trough at 70-80%.

Larger drop behind rim requires stronger B-field, in general.  Thinking about
how to formulate something quantitative out of this work.  The drop behind the
rim tells us, the damping / B-field must be at _least_ so strong.  Two effects:
thin rim (width of the bump) tells us damping is fast.
large drop behind rim tells us whether shock field is strong (larger contrast)

We should separate these two effects, which both contribute to a perception of
"sharp" thin rims.

Between 220 and 240 arcsec.
* does a max occur?
  - if no (continuous rise over 20"), damping is slow (ab large)
  - if yes, is there a trough? (damping is moderate to strong, ab small)
    + if trough, how deep is the trough?
      deeper trough means stronger field
    + if no trough, how far does emission drop?
      steeper drop means stronger field

We could repeat this for a range of values.  But, actually it's better to take
a moderate size range and see if that is enough.  Going back more will be more
and more inaccurate.

For each point in a grid of B0, `a_b`, run this simple algorithm on a grid.

Idea 2: for better X-ray modeling (if we use only 4-7 keV), can we integrate
model outputs over 4-7 keV?  Image from `merge_obs` is photon/cm^2/s, whereas
the output from the full model is (I think) erg/cm^2/s.
Then, we'd have to integrate I(E)/E^2 dE to get counts.

(probably a small effect. you considered this before)


(Week 26) Monday 2014 November 24
=================================

Summary
-------


Morning meeting
---------------

Reviewed (individual) profile fitting in X-ray and radio.  Some points
raised by Rob/Brian:

* What happens at smaller `a_b`? How about `a_b` values close to FWHM fit #s?
* What happens in the loss-limited case (even if radio rims incompatible?)
* Brian, on fit quality + B0 values estimated from indep. X-ray/radio fits:
  factors of 2-3 are nothing to lose sleep over (to fit the astronomer cliche)
* We don't know the data resolution for Steve's radio data (what configuration?
  good question)
* Joint fitting of profiles looks like the best approach
* In joint fits, check difference of best fit translation factors and ensure
  the number looks reasonable.

Chandra mtg: rooms at Boston Park Plaza are small and old, apparently
Brian talked w/ Steve, Frank Winkler, Knox Long... update to SN 1006 results is
a big question.  We'll also want profile shapes/fits for SN 1006, in whatever
comes out of next SN 1006 project work.

New SN 1006 EVLA data (13A-154, PI: Dave Green, L band, 2013/2014) pending.
Dave Green met w/ Rob last week -- EVLA is a whole new instrument, scientific
community is still getting acclimatized.  Need some time (weeks to months?) to
properly reduce and get scientifically usable data out.

Notes to self:
* compare model profiles from best _FWHM_ fits (old procedure) to data (!)
* in fit plots, show all data (included and excluded from fitting)


Radio profile shapes
--------------------

Using code (Saturday Nov 22) to calculate useful profile parameters:
* does profile have a max?
* max position?
* does profile have a trough behind rim?
* min position?
* min intensity?

Note that (does profile have trough?) is equivalent to
(min position != rminarc), so this is redundant.

Very interesting result -- if the trough behind the rim is resolved (by rminarc
setting), the (relative) intensity drop from thin rim appears to be
_independent_ of the damping length, yielding a nice curve of intensity vs. B
field!.  If the trough is not resolved (weaker damping, larger `a_b`), we see
curves splaying outwards from expected curve (smaller intensity drop at given
B0 than expected.

Why?  There is, I suspect, a self-similar-ish solution to be had.
Something that depends on B(x)

We also have nice curves of the trough location (if resolved).
At small B field, the trough is very close to the rim for all damping lengths;
as B field increases, trough move backwards.  Trough moves back faster for
weaker damping (larger `a_b`), sensibly.

Cleaned up some preliminary/prototype plots and code..


Tuesday 2014 November 25
========================

Summary
-------
* Profile shape, bounding params from minimum intensity (plots, test etc)
* Code to fit x-ray/radio profiles jointly, some tests

Profile shape analysis
----------------------

The idea: see how far emission falls behind the thin rim, and use this to set
a stringent lower bound on the magnetic field (and, a less certain upper bound
on damping length).  My assumption is that relic / old emission _falls off
monotonically_ towards the forward shock.

Put together plots, material on constraining magnetic fields / damping length
by using this new approach -- measuring the relative drop-off in downstream
emission just behind rim (steeper yields better constraint).  Ran by Brian --
just write it up, so they can all look at it and assess whether it's helpful.

With an ad hoc script (doing by hand, essentially), I can extract the
intensity minimums in X-ray and radio (large, new regions for VLA data), and
the position behind the shock, roughly.  We need a way of pulling out the
position + maximum drop (steepest drop).  Get a few data for each profile, and
compare -- which datum gives the most stringent limit?

This won't be that quantitative, because a whole subspace of parameters will
work within our intensity bound (sufficiently small damping length and large B
field will always suppress downstream intensity).  But, if the model is to
exactly (pipe dream) match the observed shock falloff, we will favor a much
smaller range of parameters (trade-off between B field and damping length).

Right now, what I'm able to do is to grid over B0 and ab, as finely as I'd
like, for an arbitrary rminarc value.  I think the procedure should look like,
for radio profiles:

1. pull out downstream limits from profiles, I suppose by "hand" for now
2. compare to contours of various `a_b` values, and find roughly allowed values
   of magnetic field and damping length.

Another thing we can do, if the emission rises back up (in the model), is to
use that to constrain things even further.  If emission must be suppressed
behind the trough/minimum then we have to throw away more of parameter space.

This needs a separate analysis/plot of its own.

Profile fitting (joint)
-----------------------

Finished code.  When running fits, need to set the relative weights
appropriately since the units are so different.  Currently, I'm using the
average intensity...

Fiddling around last week, I often had to bound `a_b` above to get the fit to
converge (instead of running away towards higher B field and higher `a_b`) --
but that depends strongly (!) on the applied profile cuts and starting
parameters.  The way to go, it seems (from Monday meeting), is to fix `a_b`.

Kind of icky but I am getting the plots out that I want.  Now how do I get more
useful numbers/results out?

Tidbit for plotting VLA data
----------------------------

Throwing here... temporarily? who knows.

    # Note: scale/cmap settings not right for 1-1.7 keV image

    ds9 -rgb \
        -red ../VLA/TYCHO_IF1.FITS \
            -scale limits 0.00015 0.0035 \
            -cmap value 2.25 0.55 \
        -green 1-1.7kev_mosaic.fits \
            -scale limits 7e-9 4e-7 \
            -cmap value 3.59 0.12 \
        -blue 4-7kev_mosaic.fits \
            -scale limits 7e-9 4e-7 \
            -cmap value 3.59 0.12 \
        -rgb lock scale yes \
        -asinh

        pan: 0:25:15.764, +64:08:25.00

Paper
-----

Generate new (preliminary?) figure for profile shape analysis text.
Paper, for this weekend.


(Week 27) Monday 2014 December 1
================================

## Morning meeting

Rob out W/Th/F-ish next week (!)
Send paper around by mid/late-week, so people can read over weekend.  Perhaps
set up telecon next week.

## Manuscript work

Applied (small) comments/fixes, and corrected synchrotron cut-off derivation in
manuscript appendix.  Cleaned and slightly restructured discussion for clarity.

Updated manuscript tables to regions-6.
* tab-fwhms.tex
* tab-spec.tex (already updated)
* tab-fits-loss.tex
* tab-fits-all.tex
* tab-fits-all-eta2one.tex

Small point of concern -- some fits have changed even though FWHMs have
ostensibly not changed (except regions 3/5/12/19).  The reason is likely that
we're using new tables (and different code resolutions).  But the results look
generally consistent.

NEED TO REGENERATE NUMBERS FOR weak/strong DAMPING regions...!
Since we're using new FWHMs -- I think I did already, but double check

dof's for excised spectrum fits are incorrect (!).  BUT, doesn't matter since
we don't show the table of downstream spectrum fits w/ lines fit or excised


Tuesday 2014 December 2
=======================

## Paper updates/writing (cont.)

Twiddle and add spectral variation figures to paper.
Finish recomputing and writing up spectral variation results.

Regenerate plots:
* Profile and spectrum fits (`prf_*.pdf` and `spec_*.pdf`)
  (should be same as Tycho regions-5, but replaced anyways)
* Intensity profiles, model fit energy-width curves, grid of curves for all
  regions, plot of `m_E` values
* Spectrum modeling, damped case (eta2=1)

Computed spectral variation for damping fits with eta2 free.  Also split code
off to new notebook (and renamed notebooks a bit more sensibly...).

    Region	rim	down	delta	ab	B0	eta2
    1	-2.494	-2.563	0.069	0.008	25.4	1.63
    2	-3.257	-3.370	0.11	0.003	16.8	19
    3	-4.586	-4.693	0.11	0.003	14.1	184
    4	-2.511	-2.581	0.07	0.004	25.2	1.72
    5	-2.227	-2.311	0.085	0.003	27.6	0.487
    6	-1.633	-1.660	0.028	0.003	114.4	0.000356
    7	-1.861	-1.914	0.053	0.004	241.6	0.0172
    8	-1.973	-2.026	0.053	0.003	52.2	0.152
    9	-1.711	-1.760	0.049	0.002	46.8	0.0176
    10	-2.151	-2.215	0.063	0.002	36.1	0.363
    11	-2.047	-2.137	0.09	0.002	31.0	0.2
    12	-1.831	-1.916	0.085	0.002	35.0	0.0518
    13	-1.771	-1.825	0.055	0.004	243.4	0.00072
    14	-1.717	-1.758	0.041	0.003	89.9	0.0119
    15	-2.814	-2.865	0.051	0.003	31.3	4.48
    16	-2.248	-2.298	0.05	0.004	340.5	0.22
    17	-2.264	-2.297	0.033	0.003	98.4	0.463
    18	-9.636	-9.766	0.13	0.006	8.7	4.1e+03
    19	-1.746	-1.802	0.056	0.002	45.8	0.0241
    20	-2.538	-2.577	0.039	0.003	47.3	1.7

    Mean 0.0660162179528
    Stdev 0.0270739322474
    Max 0.130300989869
    Min 0.0278355858292

## Radio rim analysis

I have to clean up my code / workflow for the radio data.
This has to be done in a more systematic way.
* generate new regions, better aligned w/ the X-ray selections so we can make
  more useful statements about damping / lack thereof, especially in relation
  to our FWHM measurements
* make plot of profiles for selections, in radio + hard X-ray
* identify steepest emission fall-off behind shock, if present
* do the same in X-ray and see what that gets us when plotted...
* Add more, whatever else is possible to better constrain profiles.

But, need to send manuscript around first

Merged (very preliminary) VLA data/region set-up into regular `data-tycho`
structure, roughly.  The derived data products will of course be different for
the radio/X-ray region selections.  Generating new set of regions.

Can we use the Dickel / Moffett / Hewitt VLA data to trace how structures,
knots are moving in radio data?  Clearly some synchrotron emission is
correlated w/ X-ray thermal emission, at knots and filaments.  Filaments are
understood (we think), but what about the knots?

Selecting new regions -- considerations are
1. avoiding radio "structures",
2. getting enough x-ray photons in 4-7 keV

Get proper motions from Katsuda as a function of azimuth angle and then
generate a table of 15 yr shifts (estimate..)

Note: in joint fits, there is some tension in trying to fit radio and X-ray
profiles simultaneously -- perhaps suggests that allowed parameter ranges
conflict? (if, e.g., x-ray rim is too sharp, but radio rim doesn't allow for
much damping)

Wednesday 2014 December 3
=========================

## Chat with Brian about fits and writing material

Where are we going with these fits?
Fuzzy, chi-by-eye... we're going for a plausibility argument.
Can our model match what we see?

(Brian [would have] thought we'd already chosen the appropriate regions, the
first time around... get it right the first time, or move quickly out of the
prototype stage)

It will not be a precise science (vs. pulsar timing or what have you), just get
used to it.  We shouldn't be overfitting, since we know so little about what's
going on (all the stuff along the line of sight, plus projection effects, plus
stuff behind too -- exaggerating slightly).

Find some kind of sensible fits -- however you like to do it.  Can we
reasonably make it work?  Use intuition.

Put together best possible plots.  Show cases where it seems to work, show
cases where it only works for one (radio or xray) and not the other.  As
always, don't show all...  Show what happens when we plug in reasonable
parameters.

## Radio profile analysis pipeline

Start with `regions-6-VLA`, 12 regions to test.  Generate profiles (radio and
4-7 keV), convert projections to boxes for overlay image on SNR.  Generate az
angles and interpolated proper motions as usual.

__Proper motion__: use data from Katsuda et al. (2010) -- favor over data by
Reynoso+ (1997) since our radio/X-ray data are from 1994/2009 respectively.
Alignment definitely looks better w/ correct proper motions.
Currently using nearest neighbor interp, seeing that Katsuda+(2010) averaged in
each sector to get a mean proper motion.

__Profile errors__: use RMS error from radio data; for X-ray data, use count
errors as usual.  Since this may be a different image than that of
Reynoso+(1997), the error could be larger.  The RMS I see in an annulus around
the remnant is ~0.0003.  This seems reasonable, even looking at smaller regions
outside of the remnant (sampling darker/brighter regions).  So the est.
error is ~0.3 mJy/beam, about 3x larger than 0.110 mJy/beam stated by Reynoso+.
(which makes sense, if this image is only A configuration).

## Radio profile, region selections

__Region selection__: After updating procedure to align radio/X-ray profiles,
I extend radio profiles to extend back at least 10% of shock radius (24 arcsec)
into remnant.  Threw out two regions that aren't necessary.  Adjust Region A
(SSE, not coinciding w/ our original regions-6) to not overlap w/ a radio pt
source behind shell emission (RA 00:25:26.02, dec +64:04:40.4; J2000)

__Why pick new regions?__: because we need more counts for profile analysis,
vs. just fitting for a FWHM, we would like more counts.  Also, we no longer
need to resolve the FWHM at lower energies.  There is of course a tradeoff in
taking larger regions, that we're averaging very heterogeneous behavior behind
the shock (various knots, clumps, etc in radio and xray).

__What's up with region 18?__: we don't select a radio/xray region in region
\#18 (of Tycho regions-6) because the signal is faint, and the orientation
changes a lot between radio/xray over 15 years.  Area appears somewhat special
and may be associated w/ ejecta outbreak just behind.

If the ejecta outburst is driving a new and very rapid shock into the ISM, or
has somehow re-energized the shock, and the old radio image traces the same
shock, then it's moving almost twice the speed of the rest of the remnant
(~8000 to 10000 km/s).  Not rigorous at all, but this could plausibly explain
why Region 18 shows a stronger fall-off in rim width.  Very qualitatively, the
old radio image also appears to show a weaker / no rim around that area.

__Selection bias__: we are definitely picking regions that have radio rims,
over those that do not.  Why?  Simply, we can place better constraints on those
areas, and they are typically associated with clear X-ray emission too.
But we should take regions all around the remnant, to say something
about damping/loss-limited behavior as a function of azimuth angle...

As much as possible

PLOT CONSTRAINT B FIELD AS FUNCTION OF AZ ANGLE... THAT WOULD BE USEFUL.

What is the effect of ejecta clumps on radio synchrotron emission, particle
acceleration?  Would the magnetic field be locally stronger, deformed around
the clump?

Structure near rims could give rise to spurious rims / fall-off, complicating
interpretation.  I try as best as possible to stick to filament-y things



Thursday 2014 December 4
========================

Orion EFT-1 launch scrubbed, today.

## Radio profile pipeline

Extract and shift relevant profiles.  Hand pick cuts for fitting profiles, and
locations to extract min/max emission (allow numpy to identify max emission if
possible, otherwise eyeball and select max location by hand).

Generated plots of % emission drop vs. proper motion and azimuth angle -- no
obvious correlations or anything.  Looking at total emission vs. azimuth,
there's clearly some trend, but that's sort of obvious looking at the remnant.
So the most useful analysis is probably to tie results to X-ray FWHM fits and
qualitative observations about the shock structures -- associations with ejecta
knots, comparison of X-ray/radio brightness, dynamical age of shock (suggested
by proper motions).

## Fitting profiles, again

Fitting profiles -- as I noted before, often fits don't converge nicely.
Models predict too steep a rise at the shock.  Fitting keeps running and
hitting up against translation factor limits.  I don't really understand why...

A little bit of reading on fitting / optimization problems:
* MCMC for Bayesian inference, e.g. by Kwon, Looney, Mundy (2011, ApJ),
  (http://adsabs.harvard.edu/abs/2011ApJ...741....3K)
* Ad hoc penalty likelihood functions, e.g. Bissantz and Gerhard (2002, MNRAS)
  (http://adsabs.harvard.edu/abs/2002MNRAS.330..591B)
But these are sledgehammers for a nail.  Interesting, and lots of material I
need to learn to use someday...

One possibility might be to penalize the fit if it moves the data beyond the
model computation domain (increase the safety factor though), but since we
can't resolve the slow rise / precursor at the shock we don't want that to
happen.

Another is to not let the model rise above the observed emission at any point,
since that is unphysical.

Perhaps try generating a small grid.  Then twiddle the factors by hand on the
ones that look passable.  Save the grid parameters to a file (just a JSON or
something).

I forgot that overloading the IPython notebook w/ console output will quickly
slow things to a grinding halt.  Remember to use `%%capture` for that, and save
your output to disk in case anything breaks.

Ran by Brian -- send today or tomorrow, okay.  Meeting on Monday as usual...

## Manually fitting (and muttering away)

Looking at Region A.  The best fit I'm seeing is around ab = 0.02, B field
around 50 microGauss.  But there's a lot of leeway in the fits.

`$B_0 = 75$`, even `$100$ $\mu$G` looks passable.  Good amount of leeway in
fitting rims, sort of tilting/sharpening rim edges slightly.  But, `$B_0 = 30$
$\mu$G` looks bad -- starts to need a new mechanism to generate rims (ejecta or
whatever).

I didn't fit all the profiles yet, but just by eye:
* Regions B, C, D, maybe O permit weak damping.
* Regions J, L, M permit a loss-limited field
* I think O is associated w/ the NW ejecta breakout; shock changes direction visibly over 1994–2009

Sent email with some example profiles.  Notebook littered w/ cells testing
varied parameters.

Friday 2014 December 5
======================

Out for day.  Minor twiddling of radio/X-ray processing notebook.

Saturday 2014 December 6
========================

Settled on an approach to constraining `a_b` and `B_0` based on rim drop-off
(and steepness of said drop-off).  We could do one
Set up, I think, a working approach the simplesto generate

Updated full model code error messages to be more customizable
(and made some of them shut up).  Also allowed the FWHM calculation step to be
skipped -- saves maybe 10% time.

Set up the profile shape discrimination.  Realized it must be an apples to
apples comparison (how steeply does it drop to _observed_ emission minimum, not
the _modeled_ emission minimum).

Next question: what should we do if a rim is NOT present?

Monday 2014 December 8
======================

Meeting
-------

Show the plots w/ varying B field
One paragraph on shape analysis.  State as a soft bound, language not so firm.
We don't need to show all the regions.  Fits are not to be trusted.
We cannot argue for good/bad data.

What happens if the X-ray profiles are used to constrain shape as well?
(in my maps of bounded parameter space) ... if not done don't worry about it.

Rob out Thursday, Friday (one day meeting), back Monday.
Needs airplane reading.


Analysis, pipeline, eyeball fits

Finished setting up soft "bounds" analysis for all rims, but I didn't
incorporate X-ray bounds yet.  That's a small step but given that this shape
analysis may not be useful, probably not worth going as far...
It's been two weeks already too...

A general observation -- IF the magnetic field is large, ramping the B field
up/down and then rescaling the profile does not seem to do much, especially in
radio.


Tuesday 2014 December 9
=======================

Summary
-------

* Finished eyeball fits for all `regions-6-VLA` selections
* Generate plots of multiple fits for manuscript
* Throw figures into manuscript, update text for radio/X-ray fits
* Sent new manuscript (finally!) to Rob/Brian

A ton of parameter twiddling today, and finicky details on plots... time to
finish writing and move on to other things.

Eyeball fits
------------

I reproduce the numbers here without comments to be sure, but they are stored
in `data-tycho/regions-6-VLA/fit-params.txt`.

    Region  1 (A): B0 =  50, ab = 0.02               DECENT rim example
    Region  5 (E): B0 = 120, ab = 0.03               bad, 2 flmt?
    Region  6 (F): B0 = 300, ab = 0.025              bad, 2 flmt?
    Region  7 (G): B0 = 250, ab = 0.02               DECENT
    Region  8 (H): B0 = 300, ab = 0.02               DECENT
    Region  9 (I): B0 = 250, ab = 0.02               bad radio rim
    Region 11 (K): B0 = 400, ab = 0.01               GOOD
    Region 14 (N): B0 =  23, ab = 0.004              bad, 2 flmt in NW
                   B0 = 800, ab = 0.01               alternate fit #s
    Region 15 (O): B0 = 200, ab = 0.005              bad
    Region 16 (P): B0 = 150, ab = 0.012              bad, multiple flmts?

    Region  2 (B): B0 = 200, ab = 0.05     (plateau) GOOD plateau example
    Region  3 (C): B0 =  15, ab = 0.01     (plateau) bad, few x-ray counts
    Region  4 (D): B0 = 100, ab = 0.05     (plateau) bad, same as above

    Region 10 (J): B0 = 300, ab = 0.3      (rise)    bad rise, doesn't fit well
    Region 12 (L): B0 = 250, ab = 0.2      (rise)    DECENT
    Region 13 (M): B0 = 200, ab = loss-lim (rise)    GOOD

So I think Region B (plateau), K (rim), M (rise) are the best examples.
I is a nice example of where the radio rim does NOT work.

Every region has a little tension, unsurprisingly.  But only Region N is
especially bad.  The rest are slightly off in X-ray rim shape, or radio shape
is funny/bumpy.

Overhauled fitting notebook / cells to do what I need.  Now as compact as
possible (I think), and generates manuscript ready plots w/ varying ab and B0.

I reached 1789 cell evaluations before resetting for today (over Mon/Tues), and
I think I got about 200+ more in after resetting the notebook.  Probably 2000+
cell evaluations over two days.

Manuscript
----------

New figures, new text.  Updated conclusions and abstract but it is pretty
choppy right now.

Wednesday 2014 December 10
==========================

## Advection solution check
Quick review of 1st order PDEs.  Linear non-const. coeff. PDEs are all just
some twisted form of the advection equation, so transform accordingly and
solve.  Confirmed Sean's solution (and, analysis is good for damped and
undamped cases alike).

To incorporate a self-similar velocity profile, we just incorporate it into the
function `z(x)` as:
\[
    z(x) = \frac{1}{B_0^2} \int_0^x \frac{B^2(x')}{v_d(x')} dx'
\]
recalling that `B_0` is magnetic field at shock, `v_d` is downstream velocity,
and `x` is downstream radial coordinate (i.e., sign opposite to `r`).
Then we simply have a new parameterization of the characteristic.

All verified and looks good to go!  The advection solution is very elegant.

## Paper

Reviewed paper again, making several minor edits to prose.  Throw out table of
`srcutlog` fits as suggested by Brian.

Generate diff as (after getting the 2014 Oct. 31 version of paper from git):

    latexdiff paper-tycho-spr-sr-last.tex paper-tycho.tex --math-markup=3 > paper-tycho-diffed.tex


Thursday 2014 December 11
=========================

Mid-day snow flurries!

Next steps?
-----------

Ran by Brian... spectral variation analysis is very iffy, honestly.  Probably
not worth exploring further, maybe don't even need it in the paper.  Will
discuss further next week.  Probably not worth exploring any more alleyways at
this time, we have to cut it off somewhere.

Radio data?  Don't worry about Panutti's observation (TDEM0200).  Probably not
useful.  Proper motion needs A config and they already have a 30 year baseline.

Draft AAS poster for Monday meeting would be ideal, print before leaving.

Proofread the manuscript, refine prose / fix typos etc.

Shape analysis, X-ray + radio
-----------------------------

Computed profiles for radio + X-ray shape analysis.  Similar to radio -- we get
slightly better constraints, but again all this analysis is very wishy-washy.
Can't really constrain fits, it's just one proxy for shape.
(continue this.. then address 2 muG eyeballed fits)

The plots (and code) are really ugly -- I can make nice contours (just need to
convert from lists of tuples, to meshgrids / 2-d arrays and then assign values
based on shape match or mismatch), but at this point it's probably not worth
it.  So I'll leave it as is.

Radio/x-ray fits with Bmin = 2 microGauss
-----------------------------------------

In short, there appears no obvious effect.

Region C, with ab = 0.01 and B0 = 15 muG, would have to be updated to about 6-8
muG instead.  But the other regions seem okay.  Taking B0 = 6 muG reproduces
the profile almost exactly, unsurprisingly; the ratio B0/Bmin may well
determine the shape here.  Smaller damping length will also have bigger effect,
of course.

It takes about 4 or so damping lengths (falloff about 90%) to get close to
asymptotic value.  So before then it doesn't matter much at all.

Actual data / manual "fits":

Region A: ab = 0.02, B0 = 50 muG        far downstream (behind peak) emission
                                        is visibly changed.  30 muG now fits
                                        radio rim slightly better, but the
                                        difference is pretty marginal
Region B: ab = 0.05, B0 = 200 muG       Looks almost identical.
                                        Damping length is large, magnetic field
                                        is large, so not much effect.
Region C: ab = 0.05, B0 = 60 muG        about the same
          ab = 0.01, B0 = 15 muG        change to 6 muG
          ab = 0.005, B0 = 15 muG       change to 6 muG

Region D: ab = 0.05, B0 = 50/100/200 muG    Region is already ill-constrained
                                            very marginal decrease in
                                            downstream emission at 50 muG,
                                            possibly extending plausible param
                                            range to ~30 or 40 muG
                                            But, 100 muG case looks fine.
          ab = 0.02, B0 = 15 muG            Change to 6 muG and all is fine

Region E: ab = 0.03, B0 = 120 muG           Can change to 80 muG, but worsens
                                            X-ray rim.  I really can't tell
                                            much difference.

Region F: ab = 0.025, B0 = 300 muG          Bad fit already (2 filament)
                                            w/ smaller Bmin, barely any
                                            difference.

In summary: don't need to worry about changing Bmin from 2e-6 to 5e-6
Briefly explain in text, that qualitative fits at small damping length and
magnetic field, will simply favor correspondingly smaller fields.  But compared
to the uncertainty in our chi-by-eye this is nothing.

Profiles from width-energy fits vs. data
----------------------------------------

The code I wrote for profile fitting/manipulation was really poorly designed.
Oh well.

Regenerated profiles, I was accidentally using an old version of `fglists.dat`.
These should be updated in new manuscript (but, the difference is likely
negligible).

Verdict: I'm actually surprised, the profiles and data do not disagree too
badly!  Generate new plots, just because it's a nice sanity check and the
figures are pretty.


Friday 2014 December 12
=======================

Poster making, a lot of layout fiddling/experimentation.  Continued over
weekend and settled on a nice layout on Sunday.


Monday 2014 December 15
=======================

Summary
-------
* Resend new ver. of parameter grids around (in re Steve's email)
* Final(?) poster figures + manuscript figure of radio regions


Meeting notes
-------------

### Steve's email

Counter-intuitive results = ?.  Parameter grids should cover response nicely,
send today, then send poster before Weds so Sean/Steve can review too.

The grid of model parameters might be useful for paper, to give quick/easy
intuition.

Discussed / explained weird `eta_2` effect w/ Rob/Brian.

### Miscellany

SN 1006 section -- still kind of TBD.  Can bring up?  But Brian/Rob okay with
it as is.  Setting Bmin = 2 microGauss has negligible effect (last Thurs
result).  Just put a one-liner in paper (done).

Rob/Brian out after this week.  Rob out until next year; not attending AAS
(been too hectic).  Brian back week after Christmas, then out to Seattle.

### Poster remarks

* Play with colormaps for 4-7 keV image, maybe invert (suggested by Brian).
* Several poster remarks that I implemented today.  Image scales should match;
  colormaps on 4-7 keV and 1.375 GHz images should be the same.
* Cartoon for B field is okay

Send around to all when ready, looks almost good to go

### AAS things

AAS -- Brian arriving late Sunday night.

Go to all the plenaries, that's probably the thing to do.  Except day of
poster perhaps, that would be the day to skip.

AAS probably has capped out around 3000 people (DC is one of the bigger
meetings).  Quite large (although, not compared to AGU or soc. neuroscience or
APS etc.).


Poster design/layout
--------------------

Playing with a lot of layout details.  I can't get the whitespace distribution
quite the way I like, but it's good enough I hope.

A few remarks on poster design:
Title and headers: Helvetica Neue Light
Body text: Adobe Caslon Pro

Circular images for SNRs inspired by posters of ML Wong (Caltech GPS),
especially his July 2014 Mars poster.  Two large images of SNRs (RGB X-ray
Tycho + monochrome radio Tycho) should be same size, ease comparison.  Needs
scalebar.  Twiddle the colormap some...

Spent a long time fiddling images to get things to work (and meddling with
APLpy, with no luck).


Tuesday 2014 December 16
========================

Sent poster around.  Some reading / review of papers.  Need to think about:
* Kepler in radio?
* Optical filaments?
* CR literature / derivations (not up to speed here)
* Shocks / blastwave derivations

Extra model (FWHM-energy dep.) fits
-----------------------------------

Fixed bug in model code (because `Fullefflength_port.fefflen` output was
updated, the adaptive 1st pass calculation needed to be fixed).

Ran tycho loss-limited fits with eta2 = 1 fixed.  As expected, almost no
difference.  B fields change by maybe ~1 microGauss, completely negligible.
Well, it's done for completeness's sake and doesn't take very long.

### Remnant distance and damping

Effect of distance to remnant: fits were originally run with Tycho `regions-5`
on October 3, using d = 4 kpc instead of 3 kpc.  Because `regions-5` is largely
the same as `regions-6`, I let the stated results stand.

The only question is damping, which I expect to be unchanged.  If, e.g., we fit
with `eta_2 = 1` fixed, that will change to `eta_2 = 9/16` at 4 kpc.  But
everything else will be the same since we've shown that B field doesn't change,
and ab is distance independent.  Should be the same result when we fit/match
radio profiles.  So no need to really dig this way, especially since our
conclusions don't hinge much on the width-energy dependence fits anyways.

Wednesday 2014 December 17
==========================

Hitting minor to-dos, reading.

Telecon.  Slew of material in ~90-120 minutes.
Adding notes today/tomorrow, my brain feels fried somehow..

All notes added -- will annotate/mark up each item when done.

Poster
------

Implemented all of Steve's suggestions for the poster (I'm not a big fan of
poster refs but, they would be useful if poster is distributed digitally).
Additional text + references required a good bit of layout juggling, threw out
Region 1 plots for X-ray width-energy dependence.

(Done Thursday) Generated "special edition" plots for poster (change `B_0` to
`B`, show legend for Region 16's width-energy dep. plot, use eta2 = 1).


Telecon to-dos, points of discussion
------------------------------------

Transcribed from my notes as best as possible.  I try to emphasize the
deliverables, and what we wish to learn.

All Figure/Table numbers refer to the Dec. 10, 2014 version of paper that was
sent around.  Table 7 was excised from that one but contains loss-limited fits
to width-energy dep. w/ `eta_2` fixed at srcutlog derived value.

### Major points

0. (DONE) Steve has only a few remarks on the poster, will send today.  Need to add
   some references.

1. Differentiate between damped/loss-limited behavior.  Behavior is obviously a
   continuum.  Plots on top right w/ damping lengths `a_b = 0.05` are clearly
   controlled by losses in X-ray, so we can't sensibly call that "damped".
   The literature is a bit misleading in presenting them like dichotomous
   models.

Brian: large grid of width-energy dep. plots (full page thing) seems like one
of the key figures in the paper.

Compare FWHM to `a\_b` for each region to get at this.  Better, use max(l\_ad,
l\_diff) / a\_b as the proxy.  Since l\_ad/l\_diff are just functions of the
parameter space this draws a line through the space of B\_0/a\_b, delineating
damped vs. not-damped in radio and X-ray separately.

For each region, tabulate:
- FWHM at 2 keV
- a\_b from width-energy fit (mu=1, eta2=1)
- a\_b from eyeball radio/x-ray fit
- l\_ad, l\_diff for mu=1 (functions of eta/eta\_2, B\_0, v\_d, nu) in X-ray,
  radio

NEXT, to just explore/characterize parameter space, compute values over grid of
B\_0 and a\_b values.  Obtain:
- predicted FWHM at 2 keV for the parameters considered
- max(l\_ad, l\_diff) for radio, X-ray
- plot on grid.  We'd expect to see contours of constant max(l\_ad, l\_diff).
  These should demarcate transition between damping and loss-limited behavior.

2. (DONE) The plots of m\_E as a function of energy are confusing.  What are we
   showing?  Steve as asking -- if we got the profiles at low energy then they
   should show that behavior (or, I would rather say that the model predicts
   that behavior, which would be affirmed or ruled out at low energies, 10s of
   eV being in UV now).

(DONE) Make a plot of these "weak-field damping" profiles for multiple energies.

3. (DONE) Lower minimum B field to 0 microGauss.  Do this for:
   - model profiles (giant 20 plot grids)
   - width-energy best fits (eta2=1, mu=1)
   - eyeball x-ray/radio fits

Took notes on eyeball fits, made new tables for width-energy fits, made new
parameter grid.

4. (DONE) Sphericity?  Brian asked about the equations and assumed magnetic
   field behind shock, etc.  All planar assumed.  Mention early in paper (model
   exposition) that we are neglecting these effects.

5. Discussed/explained the "effective velocity" explanation for D(x) behavior.
   Wouldn't D(x) be very large far downstream?  Could explain by electrons
   being diffused (large mean free path) so fast that they don't contribute
   meaningfully to the rim, smeared out over huge lengthscale.

Effect is still present without D(x)B^2(x) = const. assumption, just not as
strong.  Should note that this is not super-exponential because (e^(-x))^2 !=
e^(-x^2), big difference (!).

As far as the velocity effect goes, the sharp change should only be relevant
right at the shock -- downstream it should be as normal.  I think that's still
okay -- we impede motion near the shock, where things are energetic + diffusion
is weak, just reinforcing the effect.  Farther downstream, everything is back
to normal and particles move rapidly away, reinforcing the thinning effect.

One question: what would happen as B\_min -> 0?  D(x) will blow up towards
infinity, rather than a finite downstream value!

6. Veered back to discussing l\_ad / l\_diff.  See my notes on 1. above, which
   combine everything together.

7. Discussed reasonableness of extreme eta\_2 values.  Do we have a
   consistency check on the eta2 values?  Steve getting flak at COSPAR about
   eta2 < 1.  Answer: srcutlog table -- put that back for Steve and Sean who
   want to see it.

Discussion of how `E_cut`/`nu_cut` scale with `eta_2`: although they're assumed
to scale inversely, Steve mentioned Jokipii (1987) suggested that with some
effects, anisotropic diffusion etc stuff I didn't quite catch, it could scale
together.

Toss Table 5, eta2 is so jumpy that it's not helpful.  Brian asks: can we do
fits with `eta_2`=10 fixed instead?  To be closer to the srcutlog derived
values.  So do this for Table 6 (also, use `eta_2=10` or `eta2_srcutlog` for
damped fits).

8. Notation: we should call magnetic field B (or B\_1 or whatever), not B\_0.
   As typically B\_0 is taken to indicate the upstream, ambient / galactic
   magnetic field.  Inherited unfortunate notation from Rettig/Pohl whatever.

9. Spectral variation discussion.  Ambivalent feelings, I think.  Steve is
   interested in the discrepancy, why doesn't our model replicate the spectral
   softening?  Try feeding in larger `eta_2`, use the values from srcutlog
   (being physically motivated), see what comes out of the model spectra.

Steve expects that larger `eta_2` will give larger `\Delta\Gamma`.  Aaron: I
wonder if that just changes indices without affecting softening/hardening; if
anything I'd think diffusion should make spectra more similar, unless strange
dD/dX effect is at play!

Separate the sections -- the organization is confusing.  This modeling was not
related to the `srcutlog` fits.


### Other remarks:

* Sean: could bring out even more the result, that radio rims require damping!

* Timeline of work: now -> hit TODOs -> send new draft -> discuss at AAS
  (Steve, Brian, Aaron) -> Steve can go "wordsmithing".

### SN 1006

* Steve figuring out what to do (consult w/ ApJ) on Sean's paper.  The
  interpretation of our results here bears on SN 1006 as well, of course.
  Knox Long very concerned.

* Radio rims in SN 1006?  Steve has been pointing out how radio rims require
  damping for a while?  Looked at Steve's conf. proceeding (Roger and
  Landecker, 1980-something) on SN 1006, cite this on why thin rims require
  damping.Actually, see Reynolds and Gilmore (1986, AJ) which has a nice radio
  image + rim profiles (http://adsabs.2harvard.edu/abs/1986AJ.....92.1138R).
  Cite this in discussion of SN 1006, perhaps.

### Misc.

* Lock Dave Green in a room with no internet to get him to do the SN 1006
  reduction, new EVLA data
* Steve's turn to find another undergrad to keep working on this rim project


Other paper/etc stuff w/ Brian
------------------------------

### Stuff on my agenda

* Intensity comparison between radio / X-ray?  Absolute flux calibration is
  hard, probably no better than 10-20%!  Radio is... turtles all the way down.
  Flux based on previous telescope, previous telescope, all the way down to
  Crab or something, who knows.  Similar issue with Chandra and other X-ray
  telescopes.  This path is fraught with thorns.  Optical is definitely better,
  but we can't use that for our work.

  Also, radio/X-ray flux changes over 15 years!  And may not change
  concomitantly.

* Moving stuff off of ACES iMac?  Talk to IT folks.  Don't use a hard drive,
  just use a thunderbolt/firewire/whatever.

* Kepler radio image?  Bottom of the barrel in terms of priority, but Brian
  will send image if he runs across one.

### Paper comments

* (DONE) Mention why we chose Regions 1, 16 to illustrate in figures (in caption),
  explain that Region 1 is one of our bad ones.

* (DONE) Figure 4 of paper -- what does this show?  The plots of width-energy fit
  profiles w/ data are better, but not sure if that should be taking up a whole
  half page...  Maybe we can just point to the grid of plots, w/ various
  damping lengths.

(DONE) To show the discrepant behavior, make the plots for weak/strong-field
damping cases at multiple energies (see 2. above).  So that replaces the few
references to Fig. 4 invoked in the discussion.

* (DONE) Tables 3, 4, 5: cut these, accept that we'll just work in the `mu=1`,
  `eta2=1` regime.  Modify text appropriately.
  (slight change of plans -- I leave Tables 3, 4 in)

* (DONE) Definition of `B_0` somewhere?

* (DONE) Text on "Formally, we do not seek a best fit value..." but we are,
  just with very coarse resolution.


My own notes/remarks
--------------------

* (DONE) make everything (all plots) based on eta2=1 best fits, so we can point
  to tables.  Throw out loss-lim / free eta2 fit tables, merge in
  eta2=10/srcutlog tables w/ Table 6 somehow.

* (DONE) Add one of the megagrid of plots (eta2=1) to show parameter space,
  partially superseding Figure 4

* (DONE) Cite Eriksen for the Tycho obs in paper, too! (following Steve's
  poster comments)

* Emphasize difference between fit results / parameter values from 1st/2nd
  halves of paper.  What's the significance of these weak-field damping fits?
  (addressing point 2. above)

* My damping fits only consider `a_b` < 0.01 right now, but we've seen from
  radio/X-ray eyeball fits that a\_b ~ 0.01 is actually ideal for Tycho.
  Please increase that bound to allow more fits?  We might get interesting
  results and see less of this weak-field behavior.

It might be necessary to make tables w/ more `a_b` values?  We'll see...


Thursday 2014 December 18
=========================

Summary
-------
* Slew of figure/table/paper updates

Spent a nontrivial bit of time organizing all my notes from yesterday, ton of
text just on _what needs to be done_.  Below I log and discuss the various
changes/updates.

All figure/table numbers refer to 2014 Dec. 10 version of paper.

Notes/summary/etc
-----------------

Change of plans: I leave Tables 3/4 in, so we can discuss whether to explore
eta2 / mu / ab dependence etc.  A few paragraphs of discussion center on these
guys.

Question: should we give the table of eyeballed fit numbers?  So people can, in
principle, look at the models themselves and say, this is crap/okay/whatever.

Currently using "best" damped fits with `a_b < 0.01` still, but I also consider
`a_b <= 0.05` at some point and see how that compares.  If we use that, we
should regenerate all tables and figures.

Tinkered with adding damping lengthscales to radio/X-ray profile plots, but
overlaying 3 shades + all the curves/data is starting to get confusing.

### New fits / computationS

* Run loss-limited fits with eta2=1 at multiple mu values, takes ~15 minutes.
  Updated width-energy dep. plots to show eta2=1 fits with multiple mu values.

* Run new damping fits with Bmin = 0 and Bmin = 2e-6.  Using Tycho regions-6,
  damped tables from 2014 Oct. 23, with ab values: [0.5, 0.05, 0.04, 0.03,
  0.02, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002]

* New eyeball radio/x-ray fits at Bmin = 0.  Store this in the right directory.

    Region A: B0 = 100, ab = 0.03       (larger ab and B0) B0 range: 50-200
    Region B (pla): B0 = 200, ab = 0.05       Same as before
    Region C (pla): B0 = 0-50, ab = 0.05      Small ab fails now, can't do plateau
                                        But, ab=0.05/B0=60 case worked before too.
                                        Just too ill-constrained
    Region D (pla): ...
    Region E: B0 = 120, ab = 0.04       ab marginally larger, but difference is
                                        small.  B0 > 50 works.  Really iffy.
    Region F: ...
    Region G: B0 = 250, ab = 0.02       Looks the same
    Region H:
    Region I:
    Region J (rise):
    Region K: same
    Region L (rise): same                      Looks the same.
    Region M (rise):
    Region P: almost the same           Could increase damping length from 0.012
                                        to 0.014, but no one cares

  Some notes on these eyeball fits specifically:
  - Loss-limited cases: even if ab is finite, B0 is too large, so we don't see
    any noticeable difference.
  - Need larger ab for radio, because we need to replicate radio rims.  But, to
  compensate we must also raise B0 in X-ray.
  - Plateaus are harder to replicate now because we don't have troughs anymore.
  The best we can do is use large damping length and get gentle rolloff.
  - Weirdly, B0 ~ 0.01e-6 still gives a (fatter) thin rim, despite being absurdly
  small.  As B0 gets small towards zero, the rim seems to approach an
  asymptotic shape.
  I suspect it's probably approaching a radio-like insensitivity to B0.
  As B gets too small, the synchrotron cooling time ramps up
  (yes, because t\_synch ~ 1/B^2).
  Also, diffusion will be very large since D ~ 1/B, so e- distribution should
  be smeared at all energies...

### Remarks on effect of Bmin=0

Drastically changes lower left, as predicted by Sean -- we recover thin rims
for ALL parameter space, because no B field means no emission.

(if shape roll-off occurs at B0 = B-no-rim = f(ab) * 4x Bmin, then as Bmin->0,
shape roll-off boundary B-no-rim -> 0)

### Weak/strong-field damping

We can use a simple criterion to distinguish weak/strong field damping.  Does
the radio rim have a FWHM?  If so, that's weak field.  Agrees w/ previous
method (looking at `m_E` plot between 1 eV to 1 keV or whatever), though not
guaranteed to agree.

Looking at results w/ eta\_2 = 1 fixed:
* weak-field damped fits, eta2=1:    Regions [1, 2, 4, 5, 10, 12, 18]
* weak-field damped fits, eta2 free: Regions [1, 2, 3, 4, 5, 10, 11, 12, 15, 18]

Okay, this changes our discussion.  Also a reminder that the division between
these two scenarios is very arbitrary.

### Calculating advection/diffusion lengthscales

with aim of characterizing damping vs. loss-limited regimes/transition, roughly.

Caught a subtle error in the `l_diff` scaling, neglecting variable mu.
Doesn't really matter since mu=1 almost always, but should present equation in
text properly.


List of new content to integrate into manuscript
------------------------------------------------

* Plot of weak/strong field damping profiles for eta2=1, mu=1
  Better Figure 4 (shows lower energies) that can also replace Figure 7.
  More informative and intuitive than Fig. 7's plot of extrapolated `m_E`.

  This is redundant -- same concept is conveyed by large grid of parameters,
  this just focuses on one particular trend.  But, include for now, for
  discussion.

* Regenerate Figure 6 with eta2=1 fixed for both loss-lim. and damping.
* Re-added Table 7 (srcutlog spectrum fits + loss-limited fits)
* Added parameter grid with eta2=1 to show parameter space for `B_0` and
  `a_b`; made minor tweaks for manuscript version.

List of removed content
-----------------------

* Figure 4.  Shifted down to discussion, to show weak/strong-field damping
  instead (not really removed, but rather repurposed/improved)
* Figure 7.  Discussion modified to address new figure of damping profiles at
  low energies, instead.
* Table 5.  We now emphasize only fit results for all regions w/ eta2=1 fixed.

List of updated content
-----------------------

* Figure 5 (width-energy dep for Regions 1,16; multiple mu) updated to show
  eta2=1 results only.  Caption updated.
* Figure 6 (grid of width-energy dep for all regions) updated to eta2=1 only.
  Caption updated.

Width-energy dep., damping fits, with a\_b > 0.01
-------------------------------------------------

For Tycho we have been using a\_b values: [0.05, 0.04, 0.03, 0.02, 0.01,
0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002].  Clearly, the range
between 0.01 and 0.05 is not very well sampled.





Friday 2014 December 19
=======================

Summary
-------

* Assembled some stuff on Bmin=0 etc for meeting


Cleaning house
--------------

Scanned Pacholczyk Ch. 3 (up to polarization) + Ch. 4 to be sure, before
returning to library next week.

Cleaned ~/Desktop/ files lying around, threw away old fit tables (the
pre-tabulated FWHMs), and discarded all raw Chandra ObsIDs (~80 GB in all).
Discarded indiv spectra from specextract pipeline (files generated prior to
merge) in Kepler and Tycho data folders.
Discarded raw/repro'ed Chandra ObsIDs (~80 GB) and reprojected ObsIDs.
Discarded energy band images extracted from Chandra (Cas A, Kepler), they are
saved in the `data-*/` folders

NOTE: if I have to calculate new spectra at home, I will have to redownload the
ObsIDs (spectra must be extracted from each individually, then merged).
I can't remember if just the reprojected mosaics were needed... but I think you
need the full ObsID.

This now drops my folder to 2.8 GB, a very manageable size.

Meeting
-------
(last one w/Rob, Brian in person)

### To-dos/disc. etc

Brian: if [you] can't get the paper out until Tuesday, may as well wait until
past the weekend, no one will read it until then anyways.

Yes, just chuck Tables 3/4.  Fits w/ varying damping length are really similar
anyways...

Sanity check: eyeballed fits to radio/X-ray should give reasonable width-energy
dependence (same idea as comparing width-energy fit params' model profiles to
actual profiles).  Check this as well -- again strictly for sanity's sake.

### Misc.

Frank, Knox, Una at meeting -- Una coming at least one day.
Knox: tends to be skeptical (concerned about SN 1006 to say the least)

HST observations of Kepler (in works by Ravi Sankrit (SOFIA/NASA Ames), Brian,
etc.).  Optical emission is mostly H alpha, maybe some radiative shock lines
(north of Kepler, south is invisible).  Low ionization states, O III / Fe II
ish or whatever.

Exploring lengthscale effects on FWHMs
--------------------------------------

Following point 1. from Wednesday's telecon -- can we compute
advection/diffusion lengths, compare to FWHMs.  Because, in X-ray, even if we
have some damping, the X-ray rims are often so thin that they aren't really
damped!

See email to Steve/Sean et al.

MAjor observation: our preferred parameter range (a\_b ~ 0.01 and B\_0 ~ 100
microGauss) falls RIGHT where all the lengthscales are comparable!  All are
important.

Current setup

TODO: on my plots, throw on the modeled FWHMs for a given `l_diff`, `l_ad`.
I don't know if this will go anywhere though.
(I.e., plot `l_diff` vs. B0 for 2 keV, eta2=1.  Plot `l_ad` for 2 keV, eta2=1.
Then, plot FWHMs for 2 keV, eta2=1 as a function of B0.  In this way we see
exactly how `l_diff` and `l_ad` contribute to the true FWHM)



Synchrotron self-absorption
---------------------------

Stupid question: could synchrotron self-absorption be relevant to Tycho in
radio?  This would terrifically mess up our results if so.

__Idea:__ brightness temp. of source can't exceed effective temp. of electrons.
Brightness temp: T of equivalent radiating blackbody w/ flux at said freq.
Effective temp: T from kinetic E of particles
If brightness temp. > effective temp., the photons would be absorbed by the
thermal population to bring it into thermodynamic equilib with Teff = Tbright
(I don't know if SNR shocks are thermally equilibrated or not. I'd guess not).
Following NRAO (Condon/Ransom) ERA, I check this...

### T\_eff

For Tycho, consider min, max `T_eff` at B = 500e-6, 12e-6 Gauss:

    T_eff = 1.18e6 * sqrt(1.375GHz / 1Hz) / sqrt(B/1G)
          = 1.18e6 * sqrt(1.375e9 / 50e-6)
          = 2e12 K -- 1e13 K

Checking backwards, E = 3kT (relativistic) gives 0.5-2.6 GeV, so we're in the
right ballpark.

### T\_brightness - compute/check specific intensity

What's the brightness temperature of Tycho, then?  Rayleigh-Jeans approx (fine
for radio) gives us:

    T_b(nu) = I_nu c^2 / (2 k nu^2)

where `I_nu` is specific intensity.  In Tycho:
* typical rim pixel: 0.02 Jy/beam, at a thin radio rim (see my paper figures)
* Green's SNR catalog: total flux 56 Jy at 1 GHz

Let's get the flux at 1.375 GHz from the NICMOS unit calculator thing.  Convert
appropriately and use spectral index -0.58.

    INPUT FROM THE FORM IS 
    Input units = Jy             
    Output units = Jy             
    Input flux =    56.0000 Jy             
    Index of the power-law spectrum as a function of frequency =  -0.580000
    INPUT wavelength =     299792. micron
    OUTPUT wavelength =     218030. micron
    Flux= 47.     Jy

Now convert to a remnant-averaged specific intensity.  How do we get the solid
angle subtended by Tycho?  Two approximations:
* pi * (4 arcmin)^2 = 4.253e-6 sr
* 2 pi * (1 - cos(4 arcmin)) = 2 pi * 6.779e-7 = 4.253e-6 sr
(so the flat sky approx. is fine)

    47 Jy / (pi * (4 arcmin)^2) ~ 47 Jy / (4.25e-6 sr) ~ 1.106e7 Jy/sr

Does this give reasonable Tycho pixel intensity?  Convert from sr to beam.

    VLA FWHM = 45/(1.375 GHZ/1GHz) arcmin = 32.7 arcmin
          = 0.0095 radians
    VLA primary beam = pi * (0.0095)^2 / (4 * ln(2)) = 0.000102 sr

    1.106e7 Jy/sr * 1.02e-4 sr/beam = 1128 Jy/beam.

Woops, clearly primary beam __cannot be right.__

    VLA synthesized theta_HPBW (A, 1.5 GHz) = 1.3 arcsec = 6.30e-6 radians
    VLA synthesized beam = pi * (6.3e-6)^2 / (4 * ln(2)) = 4.50e-11 sr

    1.106e7 Jy/sr * 4.5e-11 sr/beam = 0.0005 Jy/beam

Looks more reasonable, within 1-2 orders of magnitude.  Now inspect our image:

    Sum: 715.798 Jy/beam for all of Tycho.
    Area: 210725 arcsec^2
    Surf brightness: 0.0034 Jy/beam /arcsec^2

But to get something meaningful we need to integrate, not just sum.
Recalling that 1 beam = 4.50e-11 sr (VLA A config, 1.5 GHz):

    715.798 Jy/beam * (0.5 arcsec)^2 = 715.798 * 5.86e-12 sr / 4.50e-11 sr
    = 716 * 0.13 = 93 Jy

CLOSE ENOUGH.  Note that 0.13 beam/px is a useful conversion too.

### T\_brightness - estimate brightness of thin rim

Recall that a thin rim pixel is ~0.02 Jy/beam.  If one pixel is a source, that
source has total flux:

    0.02 Jy/beam * (5.86e-12 sr) * 1 beam / 4.5e-11 sr = 0.0026 Jy

Or, 0.0026 Jy/sr --> 2.6e-26 erg/cm^2/s/Hz/sr
(recall 1 Jy = 1e-23 erg/cm^2/s/Hz = 1e-26 W/m^2/s/Hz)
That's the specific intensity.  Then the brightness temperature is

    T_b(nu) = I_nu c^2 / (2 k nu^2)
            = 2.6e-26 erg/cm^2/s/Hz/sr * c^2 / (2 * k * (1.375 GHz)^2)
            = 4.5e-8 K

That's obscenely faint!  Is this right?
A source needs ONLY be 4.5e-8 K to radiate 0.02 Jy/beam at radio.  Our source
is "effectively" much more energetic, yet the synchrotron radiation that comes
out is not very bright.  So synchrotron radiation is not even close to being
self-absorbed in Tycho's SNR.

Current ratio of Teff/Tb is 1e12/1e-8 ~ 1e20.  HUGE.
`I_nu ~ nu^-0.58`, so `T_b ~ nu^-2.58` whereas `T_eff ~ sqrt(nu)`.  Thus,
Teff/Tb ~ nu^-2.08.  If we drop by 10 orders of magnitude in frequency (to
oscillations at Hz wavelengths... that's crazy low, ELF radio), we will hit
self-absorption.  Thus it seems the rims simply aren't that bright in radio.


Sanity check in re DSA
----------------------

D ~ E, so X-ray electrons diffuse more effectively than radio electrons.
This hinders the DSA process -- more particles escape at higher energies
(naturally gives rise to the decreasing power law spectrum)

The diffusion lengthscale works out the same because radio electrons are longer
lived, though.  DSA is most effective for short gyroradii -- which confines the
electrons, w/ small D.

