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


Main TO-DOs
-----------

Material to discuss, Monday?
* B-damping -- including electron cut-off adds a similar shift
  (I suspect Sean's estimate of `m_E \sim -0.1` may be a bit conservative,
  perhaps -0.2 is closer.  But needs validation)
  Note that my numbers disagree w/ Sean's by about 1% (FWHMs), larger than I
  would like (mainly, error is quite a bit larger than resolution error)
* Procedure for validating B-damping: my thought/instinct is to run a series of
  fits, estimate the smallest damping lengthscale ab that produces rims w/
  chi-squared within +/-1 of the best fit w/o damping or whatever.  Compute the
  relevant F-statistic.  Check a range of numbers -- it's possible that damping
  could produce a BETTER fit, even.
* Check Kepler numbers.

In order:

1. Magnetic damping! (validating)
2. Address resolution issue (Pacholczyk + internal integrals), before moving to
   more full model fits!
3. Set radio spectral index to 0.58 (s = 2.16) (note this in rsch notes)
4. regenerate tables for regions-6 fits.
5. Check out outliers by hand (why did the fits go wonky), and check manual
   error computation speed.
6. Generate Kepler tables, too.

Grouped by individual "pipeline":

* Tycho -- Cull regions-6 out of regions-5, apply blacklist, and manually port
           over derived data products as needed... (add notes in README files).
           Generate all spin-off/variant data products needed...
           (full model errors, after checking resolution + manual approach)

           Tackle the magnetic damping matter.  Set up a port of Sean's code
           and start running calculations.  Do some trial cases and compute
           `m_E`.  Figure out how to show when magnetic damping does matter, or
           does not matter.

           Try splitting 1-1.7 and 2-3 bands into smaller pieces! (low priority)
             (run `reproject_obs` and `flux_obs`, double check that I get the
              same files that Brian has sent, then make new cuts and numbers)

           For my own verification -- run full model fits w/ capped/bkg-subbed
           FWHMs, to compare to non-corrected FWHMs.  Do chi-sqr values drop
           consistently?  Does eta2 increase from ~0 values?  Possible that
           there's nothing -- but I want to check.
           (run after regions-6 is done)

* Kepler -- for regions-3: cut 1,2,8,10,11 (on basis of poor FWHM fits)
            Get spectra for regions-3
            Generate simple model fits (... redo on regions-3)
            Check full model resolution #s (ipynb)
            Generate full model table (check w/ Brian on numbers, first)
            Generate full model fits w/ errors

* Cas A -- pick subset of regions-0 for regions-1 from FWHMs/profiles/counts,
             ensuring that all regions can get decent FWHMs.
           generate new FWHM fits (1-9 keV) and get spectra of rims
           select specific regions outside rims and get spectra, need to figure
             out how much thermal contamination

* Tycho (radio) -- later, meddle with Sean's model and see if we can't
                   reproduce the radio rims.

* Tycho (2003) -- try reprocessing ObsID 3837 on another computer.  If that
                  fails email the CXC or something

Lower priority / various checks
-------------------------------

* RESOLUTION -- check Pacholczyk tables (!)
                check some internal numerical integrals
                check accuracy of Pacholczyk's Bessel function values
                (likely a tiny error -- but to be sure.  Generate
                 1. table w/ more entries from Pacholczyk
                 2. table w/ same entries, but more recent values
                 3. table w/ more entries, more recent values)
                More e- table values will slow down computation for sure...
                (ask Sean: any reason for table selections)

* Manual error computation -- test it out on one/two regions, estimate how much
  time it takes (and extrapolate).  Do this after the resolution checks.

* SN 1006 -- re-run pipeline at very end and ensure all is still consistent
             re-run pipeline on _individual_ SN 1006 regions.
             Does the sub-Bohm diffusion result still stand at that level?

* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Check all constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran.  Go back to Sean's transport eq'n and rederive.
* Check transport equation for pure advection case

* When paper is nearing finish -- re-run entire pipeline on SN 1006 as a sanity
  check.  But, not now.

* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me... (this requires the reprojected files)

### Conceptual/background/physics questions (look up and/or ask)

* When converting from particle energy E to radiated frequency nu, doesn't the
  formula `$\nu_m = c_m E^2 B$` implicitly invoke the delta-function approx for
  synchrotron radiation?  I have no idea how you'd do it otherwise though.
  The diffusion coefficient is set by a singular electron energy, which gives
  some synchrotron radiation spectrum.  Perhaps, for an individual electron
  this is not a bad approximation?  I don't understand the physics of the
  derivation of synchrotron stuff (MUST REVIEW ALL OF THIS...).
  Perhaps it corresponds to some median, mean energy, or power, or something.


Extra
-----
i.e., things that no one cares about

Write/find short scripts / hooks to git, to clear out ipynb output cells or
something... 



