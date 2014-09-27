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


Monday progress/material to review
----------------------------------

Kepler: new region selections, profile fit notebook, model fit results (tables)
        if possible -- make regions-3, then get spectra for subset of regions

Cas A: show region selections and profiles (plan: defer spectra/fits until we
       figure out how to tackle all the blobs)

Tycho: print copy of (latest) paper for Monday


Main TO-DOs
-----------

### Paper

Review and fill in results, discussion.
Double check and reference all numbers for Tycho (the radio spectral index, in
particular, needs to be explained -- don't spend more than ~1 hr looking that
up though.

### Data processing / analysis

Waiting on the full model error calculations for Tycho to finish, it will
probably take another day...

Kepler -- get all numbers and citations prepared.

* Kepler -- prepare model fits for Monday, showing all regions for now
            for regions-3: cut 1,2,8,10,11 (on basis of poor FWHM fits)
            Get spectra for regions-3
            Generate simple model fits
            Generate full model table (check w/ Brian on numbers, first)
            Generate full model fits w/ errors

* Cas A -- generate FWHM fits (1-7 keV) and get spectra, if there is time;
           run profile plots + counts by Rob/Brian on Monday.

* Tycho -- running full model fits WITH errors, to give idea.
           likely we will run it one more time, when I manually blacklist
           regions etc...

           Error finding code is awful.  It's definitely getting stuck and
           doing some stupid things.  Look at stdout, check why the adaptive
           stepping doesn't seem to be working.
           Try manual error fitting, and see how long it takes (extrapolate
           accordingly).  If it works... do EVERYTHING ELSE FIRST.  Then come
           back and compute some errors.

           Run full specextract on all bkg-2 files and re-link spectrum files.
           Look at stdout and see whether we can cut down error computation time)

           Try splitting 1-1.7 and 2-3 bands into smaller pieces! (low priority)

           For my own verification -- run full model fits w/ capped/bkg-subbed
           FWHMs, to compare to non-corrected FWHMs.  Do chi-sqr values drop
           consistently?  Does eta2 increase from ~0 values?  Possible that
           there's nothing -- but I want to check.

* CHECK PACHOLCZYK IN RE RESOLUTION -- avoid wasting computation time, if the
  resolution is too poor.

* SN 1006 -- re-run pipeline at very end and ensure all is still consistent
             re-run pipeline on _individual_ SN 1006 regions.
             Does the sub-Bohm diffusion result still stand at that level?

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

* How do you rule out magnetic damping?  E.g., looking at figure 4 with ab =
  5 percent of shock radius, seems like it would give a decent energy dropoff.
  How did Sean get that mE must be of order -0.1 for damping? (Brian: maybe
  Sean explored parameters, maybe from old/"classic" papers on damping?)


### Extra

Write/find short scripts / hooks to git, to clear out ipynb output cells or
something... 



