snr-research agenda
===================
Aaron Tran
Summer 2014

N.B. this is not an archival document, but rather replaces post-it notes and
other ephemeral pieces of paper.

LEARNING (2nd priority for now...)
----------------------------------

Get an intro level astrophysics textbook
1. basic radiative processes (thermal, nonthermal emission)
2. supernova physics (temperatures, radiation/spectra during evolution)
3. physics of Tycho -- interpretation of the fluff.  ISM gradient.
4. X-ray telescopes.  Resolution, operation, limitations

About thermal emission lines:
* [atomdb](http://www.atomdb.org/Physics/units.php)
* [snr spectra](http://www.phy.duke.edu/~kolena/snrspectra.html)
* [x-ray lines](http://www.phy.duke.edu/~kolena/strongxlines.html)

Q: why is the amplification higher than expected -- what sets the expected amt
of amplification?  (ressler mentioned expected value of 4x for unmodified,
strong shock -- where did that come from?)

(Rankine-Hugoniot conditions for strong adiabatic shock?.. go read Shu)
(Drury mentions compression ratio @ shock in DSA...)

turbulent B field amplification (nonlinear plasma physics):
doi:10.1038/nphys2978 (in a lab! this is really cool)
[ppt](http://fermi.gsfc.nasa.gov/science/mtgs/symposia/2007/p4/P4.1_Ellison.pdf)
[AR Bell](http://mnras.oxfordjournals.org/content/353/2/550.full.pdf) (highly
cited paper on this stuff)

Possibly useful link on most recent galactic SNR
[here](http://chandra.harvard.edu/photo/2008/g19/media/)

The arc in the SE -- see doi:10.1088/0004-637X/732/1/11 (evidence for
progenitor of the type Ia SN...)

Organize literature collection/review eventually...


THING in TYCHO
--------------

Try breaking it apart into different pieces, see if it looks different.

About this: see Hwang and Laming (2012, ApJ 746(2)), 3rd paragraph of
introduction:

> The X-ray-emitting Si ejecta show a bipolar structure with jet-like features
> (Hwang et al. 2004; Vink 2004; Laming et al. 2006) similar to that seen in
> optical (Fesen 2001 and references therein) and infrared emission (Hines et
> al. 2004).



Main agenda
===========

Pipeline, profile and FWHM calculation and plots
------------------------------------------------

So need to
1. overhaul specextract pipeline (handle region coordinates correctly, test and
   ensure it works correctly with my adjustments)
2. rehaul data structuring for FWHMs because it's a pain in the ass to plot /
   organize
3. include switches to test different measurements/procedures, for
   reproducibility.
4. proper config files / cmd line arguments for fitting/FWHM analysis scripts?
5. clean up rsch notes document... keep track of what changes/etc I've been
   making

* overhaul pipeline for profile fit / FWHM processing again -- it's just so
  messy, and feels hard to work with

* Ah dammit.  I think I've been using inconsistent coordinate systems etc all
  this time... so my background links, spectra, regions, etc may all be
  slightly off.  Argh. Need to fix this.

* Update pipeline documentation ...

* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.
* Shift region numbering to be more logical (instead of 1, 10-13).

* Eventually: run whole pipeline on one set of ALL regions, sampled all around
  SNR, save the output.  Use this to argue/show why regions are good/bad.

On using different calculations of the FWHM: how do I show the effect of these
different procedures?  Some kind of normalization? (Figure 10 of Ressler).
Quantify effects on calculation of `m_E`, B0, eta2.


Models for filament widths
--------------------------

* Clean up your notes.

Some high level questions / thoughts from poster session (July 31):
* Ori -- sanity check on checking luminosity and shock kinetic energy, though
  Rob/Brian note that the synchrotron radiation is very inefficient (esp.
  compared to when the remnant cools, H starts recombining and emitting thermal
  lines, T ~ 10^4 K)
* Terri/Amy -- why are the shocks so thin, anyways?  don't have a good answer
  for that...
* Our analysis, if anything, confirms that the shocks are thin to begin with...
  but doesn't say "why".  Why shouldn't accelerated particles travel farther
  back?  Function of time / remnant evolution?
* Madura -- 3D reconstruction?
* We really need to quantify the magnetic damping issue -- can we put a lower
  bound on the relevant damping scale length, given that it basically can't
  explain the rim dropoff with energy?


### General (higher level to-dos)

* Tables 7, 8 reproduced for SN 1006, then Tycho.
  As we have wanted for the last SEVERAL weeks...
  Report chi-squared values so we can compare different
  values of mu... and report errors from brute-force chi-square
  Also, print out `m_E` values, point to point and from `width_dump` model...
  see some of the old notebooks in `code-profiles/` (merge functionality
  together and throw out old / unused things)

  Plan -- use ipython notebook to prettyprint output, but use `models_all_exec`
  to abstract away as many details as possible.

  Idea -- avoid making ipython code "general".  Use it for one off scripts,
  e.g. make a new notebook for each supernova remnant.

  To do this...
  1. need a quick rootfinder for simple fit errors. (DONE... barring kinks)
  2. need a quick rootfinder for full grid fit errors.

  SOOOOO how to do the fits?

  Find the best B0 and eta2 naively in the imperfect grid.  Fit at that point,
  let it run around (~20 calls), then compute error.

  So, 6 mu values * 13 regions = 78 fits... * 20 calls = 1560 calls, 5 seconds
  each requires 2.2 hours... blahhhh.
  
  To compute error, check the eta2 grid until you find chisquared values
  bracketing `chisqr_min + 2.7`.  Perform actual fits to B0. If values no
  longer bracket `chisqr_min+2.7`, move onwards...

  (problem -- what to do if multiple values bracket... given resolution
  issues, etc)

  Tricky part is that this has to be somewhat done by hand, in case that no
  such error exists (e.g. if eta2 runs to 0 or infinity without blowing up
  chisqr...)  I don't know the best way to do this, trying to automate the
  process but still leaving room to do some processing by hand...


  If I really wanted to nail down fits...
  for each mu, for each filament, do the following:
  1. select rough region of best eta2 values, and for each eta2 run a fit
     to get the best possible B0 + chisqr value. (~20 * 20 function calls)
  2. with best eta2 + B0, run about 6 fits with different epsfcn values,
     keeping the very best fit values you can find (~6 * 20 function calls)
  3. use grid to find brackets on errors.  At each putative error point, you
     need to perform a fit to verify that that is the best error bound.  If the
     chisqr drops below threshold, move on.  (~10 * 20 function calls)

  result: 6 x 13 x (20+6+10) x 20 x 5 seconds = 78 hours total.
  Well, if I ignore step 1 and just run fits around grid best value, that cuts
  the time in half to ~40 hours.

  For SN 1006 that's 8 hours = 6 x 5 x (6+10) x 20 x 3 sec.

  This looks about right, my time estimate... annoyingly slow, just to get
  error estimates!!...  Maybe fit a curve and narrow estimates that way...
  I have no idea, urgh.
  Or, could I do it by hand?
  ORRRRR maybe given that the core issue here is that, B0 is not well resolved
  (surprisingly or unsurprisingly), I should do that parabola fitting thing to
  get a better B0 value...  That requires 1 function call, but 1 function call
  for every eta2 (for every mu).  So 1 * 100 * 6 = 600 calls * 3 sec = 0.5 hr.
  But, that's a relatively low cost... esp. if I can show it does improve
  fits... definitely faster than this error bound pushing!!!!

  Other approach -- with a priori knowledge, follow sean's approach of
  manual fitting to estimate errors (calculate chisqr on the fly), as with the
  lmfit it definitely takes some time... but that's still probably faster than
  anything I could do!...

  Could split by 3 with ad hoc parallelization (~13 hrs Tycho, ~3 hrs SN 1006)
  For Tycho, I'd have to make some of these function calls anyways, just to
  deal with the variable shock velocity gridding problem.

  Maybe we should consider looking into speeding up Sean's code


  Right now, just fitting -- 6x5x20x3 = 0.5 hours...

* ipython notebook with results of varying compression ratio, shock speed,
  remnant distance, any other twiddleables (including/excluding energy cutoff).
  Show for both SN 1006 and Tycho.
  Use this for discussion/commentary on twiddle-able knobs.

* Remember brian's suggestion (from Friday july 25): how does mE depend on
  energy? what happens if you fit a straight power law to that???

* Write code to compute azimuth angle of regions, and estimate shock speed for
  each region.

* Check all your constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran

* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Look at azimuthal dependence of B field, robustness of numbers from approx
  python model. Brian asked, does B field scale with stronger energy
  dependence? (NOT ADDRESSED as of July 30)

### General (high-level questions)

* How does Sean define unobtainable, in Table 8? E.g., for mu = 0.5 I can
  manually fit and get a chi-squared value of 6.2 (compared to ~4 for the higher
  values of mu).  Brian: yeah, run this by Sean.

* Ask Sean if there was any reason for 1sigma error in Table 7?  I think I was
  seeing larger errors from brute-force chi-square in lmfit, even for 1 sigma
  (so it seemed like the errors in Table 7 came from sqrt of diag elements of
  covariance matrix, which would be inaccurate here)

* Sean used energy cutoff in all his model fits?

### General (lower level to-dos)

* Possibly most useful: store FWHM, KEV data in SNR objects.  Maybe even this
  indices thing, for comparison with Sean's data, or to throw out 0.7-1 keV
  data in Tycho.  Call profile-fitting functions, or load pickled data, with a
  function call to initialize FWHM data...

* SAFETY FEATURE, do NOT let user call maketab if `fullmodel.so` is older than
  `FullEfflength_mod.f`... or just force it to be recompiled each time, but you
  have to reload the module in all scripts or something...

* (feature creep) write Makefiles to automate the human aspect of running this
  pipeline?

* Should we be constraining/checking resolution somehow -- checking that grid
  is "sufficiently" resolved to get a good fit?  What is the "error" of
  fractionally changing the grid... thinking of Math 228B here.

### Small checks, constants, verification

* CREATE a test case for lmfit, to verify it is doing a least squares fit as I
  would expect (just check against scipy curve fit).
  During this test -- verify that when I freeze a parameter, chi2red is
  calculated correctly with one less DOF (if it spits out a chi2red)

* Check transport equation for pure advection case
* Check numerical prefactor 8.3 TeV for electron cutoff energy
* Update numbers from Pacholczyk (?), consider adding more entries (can ask
  Sean about this)

* Check FORTRAN code for memory-efficient array indexing?
* Consider caching / memoizing tabulated electron distributions?
  Problem: advective/diffusive lengthscales depend on B0, mu, eta, etc.
  So, as we try new parameters, the electron distributions change...
  which does physically make sense, after all.  sigh.
  How can we get around the speed issue?


### Conceptual/background/physics questions (look up and/or ask)
* When converting from particle energy E to radiated frequency nu, doesn't the
  formula `$\nu_m = c_m E^2 B$` implicitly invoke the delta-function approx for
  synchrotron radiation?  I have no idea how you'd do it otherwise though.
  The diffusion coefficient is set by a singular electron energy, which gives
  some synchrotron radiation spectrum.  Perhaps, for an individual electron
  this is not a bad approximation?  I don't understand the physics of the
  derivation of synchrotron stuff (MUST REVIEW ALL OF THIS...).
  Perhaps it corresponds to some median, mean energy, or power, or something.
* How to get such strong B field values?!  Is nonlinear MHD turbulence enough?
* How far behind the shock does magnetic turbulence operate?  The equations we
  use assume electrons are injected, with power law distrib. + cutoff, right at
  the shock (I suppose as they advect/diffuse the power law distrib. may be
  modified but still persists?).  How long should the turbulence persist,
  time/lengthscale? (same, related by plasma advection speed)
* How to get larger compression ratios?  What does it imply?
* Why is cutoff energy set by equating loss and acceleration times?
  The electrons radiate faster than they can be accelerated?  But that seems
  weird -- because I thought the acceleration is what gives rise to the
  radiation.  Please read parizot et al (2006), take notes, rederive things.

* What are the qualitative differences between the models, especially in
  results?  We can see the different assumptions in source/sink terms, but how
  does that change what comes out?  We already see that the complex model
  favors smaller B fields to explain observed widths.  How does that reflect
  the underlying physics?  Complex model generates "thinner" filaments
  naturally?
* Is mu restricted to fall within `[0, 2]`, as suggested by equations for
  turbulent spectrum and diffusion coefficient?  Seems legit...
* How do you rule out magnetic damping?  E.g., looking at figure 4 with ab =
  5 percent of shock radius, seems like it would give a decent energy dropoff.
  How did Sean get that mE must be of order -0.1 for damping? (Brian: maybe
  Sean explored parameters, maybe from old/"classic" papers on damping?)

* what are the main lines of evidence for proton acceleration? Stuff Jack and
  Zeeve are working on, I suppose.
* Is there a way to relate diffusion /coefficient/ to turbulent energy?
  (wondering, what would happen if shock did not induce turbulence but only
  compressed magnetic field -- different effects for plane parallel/perp
  field, but what would those effects be?)


More supernova remnants
-----------------------

Read about work on Kepler / Cas A sometime, when you have time

* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me...
* Question: looks like `reproject_obs` will give evt file (which can then be
  partitioned by `dmcopy`, is that good enough....

Read `merge_obs` documentation, play with outputs
