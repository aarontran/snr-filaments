snr-research agenda
===================
Aaron Tran
Summer 2014

N.B. this is not an archival document, but rather replaces post-it notes and
other ephemeral pieces of paper.

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


Main TO-DOs
-----------

For discussion: can we change the lower energy bands -- or, even change them
all and aim to get proportional counts?  Look at my data on peak counts...
heavy weighting towards middle bands.

IF we add some count criterion -- make a note in the manuscript (each energy
band has at least ??? counts (integrated) at the peak's maximum datum).


* Go through `regions-all.reg` and add regions w/ bad 0.7-1 keV FWHMs.
  Review larger region selection, w/ new criteria for region selections.
  And, run them all through the fitting pipeline to help w/ region rejection.
  (regions-6, after looking at regions-5)

  IF we use the larger regions, must check the spectra carefully.

  -- Rerun the specextract pipeline -- then cull regions one more time.
  Then run the specextract pipeline, generate new plots of profiles + spectra,
  and run all fits + errors.

* Kepler -- cull regions by w/e criteria, prepare regions-2 w/ 4 energy bands.
  Go off of the peak cts/region, as I am doing w/ Tycho now.
* Cas A -- select first round of regions, get FWHMs.
  Check how we are doing on counts, may need to stick to higher bands.
  Run w/ all three high energy bands.  Might even try up to 9 keV or so.

* CHECK PACHOLCZYK IN RE RESOLUTION -- avoid wasting computation time, if the
  resolution is too poor.  (but, work on paper first... need feedback from Sean
  and Steve)


IF and ONLY IF data are worthwhile, then deal with 1. spectra, 2. full model
errors.  For the time being, only compute FWHMs, simple models, and full models
with stderr.  That's it.  No more.

EXTRA: re-run pipeline, on SN 1006 individual regions. Does the sub-Bohm
diffusion result still stand?????

### Science

From talking w/ Brian -- we're gonna need some way to show the global energy
drop-off!  Which is why he asked about any global `m_E` dependence.

* Maybe report `m_E` as estimated from MODEL fits, which, to be
  clear, is a separate "observable" that I expect to be more robust
  than the point-to-point measurements)
* Remember brian's suggestion (from Friday july 25): how does mE depend on
  energy? what happens if you fit a straight power law to that???

### Various checks

* Full model code -- check resolution error due to Pacholczyk tables (!)
  Might also be worth double checking some of the internal numerical integrals.
  And, check how accurate/correct Pacholczyk's Bessel function values are now.
  (ASK SEAN -- ANY REASON FOR TABLE SELECTIONS?)
* Remnant radio spectral indices?!
  See [ppt slides](http://www.astro.le.ac.uk/~cbp1/cta/Talks/TonyBell.pdf),
  just pick some numbers and report in text (saave in SNR catalog, w/ sources
  and explanation for choices).
* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Check all constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran.  Go back to Sean's transport eq'n and rederive.
* Check transport equation for pure advection case
* Check numerical prefactor 8.3 TeV for electron cutoff energy (DONE, YUCK)

* When paper is nearing finish -- re-run entire pipeline on SN 1006 as a sanity
  check.  But, not now.


### CIAO stuff

* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me...


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
* _How_ to get larger compression ratios?  What does it imply?
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
  (need to read more on e- diffusion)
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


### Extra

Write/find short scripts / hooks to git, to clear out ipynb output cells or
something... 



