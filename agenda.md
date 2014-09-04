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

Learning (practical stuff)
--------------------------

Maybe useful...
[GSL manual, section on nonlinear fitting](https://www.gnu.org/software/gsl/
manual/html_node/Nonlinear-Least_002dSquares-Fitting.html)


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

High level questions to consider
--------------------------------

From poster session (July 31):
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


Pipeline, profile and FWHM calculation and plots
------------------------------------------------

(on hiatus while I'm getting full model code to work and spew good results)

### Software engineering blargh

1. overhaul specextract pipeline (handle region coordinates correctly, test and
   ensure it works correctly with my adjustments)
2. rehaul data structuring for FWHMs because it's a pain to plot / organize
3. include switches to test different measurements/procedures, for
   reproducibility.
4. proper config files / cmd line arguments for fitting/FWHM analysis scripts?
5. clean up rsch notes document... keep track of what changes/etc I've been
   making

* overhaul pipeline for profile fit / FWHM processing again -- it's just so
  messy, and feels hard to work with

* I think I've been using inconsistent coordinate systems etc... background
  links, spectra, regions, etc may all be slightly off (physical coords may
  differ w/ single obsID vs. merged obsID...).  Can't remember the exact issue,
  need to double check.

* Modularize code (separate out components for profile fitting, spectra
  extraction, region processing, whatever).  Make it easier to keep track of.

* Update pipeline documentation ...

### Actual science

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

### MAJOR CORRECTION!

When performing model fits to rim widths, width errors must be 1-sigma for
correct covariance matrix scaling.

ADMITTEDLY -- this is kind of meaningless in both cases because our fits are
very often not accurate in a least-squares sense, and the errors may well be
meaningless.

For now, I DON'T apply the redchi scaling to our covariance matrix, because our
errors have meaningful magnitude; they are not merely weights.


### CURRENTLY WORKING TO-DOs

Agenda:
* debug index error (ughsdf)
* would also be nice to figure out, wtf happened to all of our logging output
* why did it disappear, are our results still okay/valid?  did it get dumped
  around piecemeal (check this by going through logs)
* start SN 1006 running, should be much quicker
* start putting together some text (set aside formatting to later)
* cobble together some kind of tycho results (decide how to store/update data)

* in SN 1006, full fits from grid + stderrs -- the tabulated stderrs appeared
  off (from a notebook run last week).  Why?!  Double check this, briefly.

* Writing - spend some time explaining why we selected what we did.
  Outline + methods within next week

* Possibly -- put out format in code-friendly format, so I can just aggregate
  all the tables together in iPython and customize the display or something
  like that?  I'm not sure.
  Could just pickle it.  I need to easily regenerate tables/plots on the fly,
  without having to redo the whole calculation again (given that it's so
  time-consuming).


passing idea: store pkl w/ best fits + chisqrs found to date...
or just maintain my own table?!  this is obscenely ad hoc.
No, too much inconsistency

Think -- how am I best going to keep track of fit/error numbers?
It is very inefficient to be constantly regenerating them (hence this mess,
and the very slow debuggin over the last ~2 weeks).  Would be good to have
stored numbers for full model unit tests, now, too.
(write some unit tests for SN 1006 -- this is just part of our validation
process)



* Perform calculations with eta2 fixed to 1 (or 0.1--10... same thing, just
  shut up and do it)

* Ask Keith Arnaud about error statistics... but, lower priority

* Generate outputs for several different data-sets.
  1. all 13 regions.  try averaging best fit parameters within filaments,
     explore the spread of numbers, whatever...
  2. arithmetic means of FWHMs in the filaments (5 data sets)
  3. geometric means of FWHMs in the filaments (5 data sets)

* Idea: maybe useful to compute advection/diffusion lengths from fits?
  ALSO we need `m_E` from this, and the data...
  (Rob, Brian: report this in the results, and make it clear
   what's going on...
   I also want to report `m_E` as estimated from MODEL fits, which to be
   clear is a separate "observable" that I expect to be more robust
   than the point-to-point measurements)

* Full model code -- check resolution error due to Pacholczyk tables (!)
  Might also be worth double checking some of the internal numerical integrals.
  But, lower priority...


### General (higher level to-dos)

* Tables 7, 8 reproduced for SN 1006, then Tycho.  This is THE high level goal
  right now.

  For SN 1006, we need 2 and 3 band fits, plots, chisqr values.
  (all infrastructure is in place, essentially)
  (most of it... just, ERRORS!...)

  Also, print out `m_E` values, point to point and from `width_dump` model...
  see some of the old notebooks in `code-profiles/` (merge functionality
  together and throw out old / unused things)

  Alternate: with a priori knowledge, manually fit to estimate errors
  (calculate chisqr interactively).

  Agenda
  5. add methods to quantify effect of varying fwhm measurements, just to
     see...
  6. Generate tables of results, for averaged filaments (using both arithmetic
     and geometric means...)

* Remember brian's suggestion (from Friday july 25): how does mE depend on
  energy? what happens if you fit a straight power law to that???

* Look at azimuthal dependence of B field, robustness of numbers from approx
  python model. Brian asked, does B field scale with stronger energy
  dependence? (NOT ADDRESSED as of July 30)

* Check all your constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran.  Go back to Sean's transport eq'n and rederive.

* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)



### General (high-level questions)

* How does Sean define unobtainable, in Table 8? E.g., for mu = 0.5 I can
  manually fit and get a chi-squared value of 6.2 (compared to ~4 for the higher
  values of mu).  Brian: yeah, run this by Sean.
  (partially answered)

  Observation:
  Sean's unobtainable values -- with two bands, his chi-squared values are
  mostly under 1 except where the fit were unobtainable.
  chi squared values range from ~2 to 16, with <1 degree of freedom...

* Ask Sean if there was any reason for 1sigma error in Table 7?  I think I was
  seeing larger errors from brute-force chi-square in lmfit, even for 1 sigma
  (so it seemed like the errors in Table 7 came from sqrt of diag elements of
  covariance matrix, which would be inaccurate here)

* Sean used energy cutoff in all his model fits?

* At what chi-sqr does our simple confidence interval analysis break down/fail?
  I don't understand the theory behind this very well.


### Small checks, constants, verification

* Clean notes, organization, code, etc for clarity
* Write some text / code pipeline explanation, cleanup, docs, modularity, blah
  blah.  Some parts are pretty good (`models_all.py`), some parts not (ipython
  notebooks, `models_all_exec.py`).
* Write some text abt methods / results so far (!!!!!)

* Check transport equation for pure advection case
* Check numerical prefactor 8.3 TeV for electron cutoff energy
* Update numbers from Pacholczyk (?), consider adding more entries (can ask
  Sean about this)


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
