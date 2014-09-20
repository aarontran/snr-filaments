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


Pipeline, profile and FWHM calculation and plots
------------------------------------------------

Tycho-specific
* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.
* Shift region numbering to be more logical (instead of 1, 10-13).
* Select slightly better regions: regions 2,3, 6,7 should extend a little bit
  more ahead of the shock for consistency.
* Slice smaller regions

On using different calculations of the FWHM: how do I show the effect of these
different procedures?  Some kind of normalization? (Figure 10 of Ressler).
Quantify effects on calculation of `m_E`, B0, eta2.

Kepler -- go through and select regions, making notes, and generate 2nd set.
And, look at higher energy bands (2-4, 4-7 keV).


Models for filament widths
--------------------------

### Main / high-level TO-DOs

Manuscript (highest priority):

* WRITE UP the transport model text,
  CLEAN UP fitting procedure text to match
  WRITE UP results section

Data, numbers:

* RUN model fitting pipeline on capped/bkg-subtracted FWHMs.
  Don't bother with full errors -- for now, just want numbers + chi-sqr values.
* PICK Tycho regions-5, according to some criteria (look at counts, spectra,
  profiles to decide).  Save regions and run pipeline on these guys.

* Kepler -- cull regions by w/e criteria, prepare regions-2 w/ 4 energy bands.
  Have simple/full model fits for Monday.
* Cas A -- select first round of regions, get FWHMs -- would be useful to
  discuss on Monday.

IF and ONLY IF data are worthwhile, then deal with 1. spectra, 2. full model
errors.  For the time being, only compute FWHMs, simple models, and full models
with stderr.  That's it.  No more.

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
* Check numerical prefactor 8.3 TeV for electron cutoff energy

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



