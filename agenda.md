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

(just to get a sense of energy -- how do we get physical units from CIAO
processed mosaics?)

Main TO-DOs
-----------

* Manuscript work (!)
  - finish up discussion.  Generate whatever new fits I need to get things to
    match up (Regions 6)

* check expressions for cut-off energy, diffusion terms in Sean's solution
  (using Bmin vs. B0 in certain places).

* Manually look at the "bad"/"spurious" energy dependence regions -- can we
  still find good damping fits that don't show this spuriousness?  I know it
  will be subjective, but just to check: if we attempt to rule certain
  scenarios / parts of parameter space out, can we still find okay damping fits
  for moderate energy dependence?

After moving forward with paper -- get srcutlog diffusion coefficients and hold
them FIXED.  Run fits to X-ray data with damping and see how the B-fields come
out for a range of damping lengths, and get radio profiles at those varied
damping lengths.

### Radio thingies

* EVLA data reduction (casaguides.nrao.edu) (prepare questions for Jack)

* Figure out how best to fit profiles in their entirety.
  One idea -- use FWHMs to line up profiles, roughly (since we don't precisely
  know the shock location + precursor may be seen ahead of shock...)
  (consider the work of Warren 2005?)
  (morlino/caprioli 2012, cassam-chenai 2007 have done this to some extent)
  Relevant links:
  http://mathematica.stackexchange.com/questions/866/simultaneously-fitting-multiple-datasets
  http://forums.wolfram.com/mathgroup/archive/2011/Sep/msg00555.html

* Radio rims -- how to get intensities to work out right?
  (Chandra photometry?... VLA already has the calibration routine)
  Can we show/prove that _thin_ radio rims are not possible in a loss-limited
  model?
  Diffusion in radio?
  time-dependence, aging, other effects may be important!

### Data stuff (just checking/solidifying results)
* Run fits to SN 1006, individual filaments (pending response/remarks from
  Steve/Sean too, perhaps)
* Set-up loss-limited Tycho fit tables w/ new default resolutions.
  Generate damping tables down to ab = 0.003, ab = 0.002, to see what happens.
  We'll be limited by the FWHM "disappearance" matter...

### Theory stuff
* The diffusion coefficient timescale depends on both upstream/downstream
  coefficients, and requires a simple integral.  What happens if we allow the
  diffusion coefficients to vary?  Does this timescale change significantly?
* Rederive the advective (transport) solution in damped and non-damped cases...


### Lower priority:

1. Use srcutlog (from Brian) to run fits with eta2 fixed at values determined
   from synchrotron break frequency.
   First make short pyxspec script to get srcutlog fits, then compute the
   desired diffusion coefficients.
2. Check outliers by hand (why did the fits go wonky), and check manual
   error computation speed. (barely relevant anymore)
3. regenerate tables for regions-6 fits (also damping tables)
4. Generate Kepler tables and run full model fits
5. ODE solving, simple model with damped B-field.  Not likely to do much, esp.
   if our full model comes up correct...


Data processing to-dos, by individual "pipeline":
(i.e., getting FWHMs/regions/etc)
(set aside while we address magnetic damping etc...)

* Tycho -- Cull regions-6 out of regions-5, apply blacklist, and manually port
           over derived data products as needed... (add notes in README files).
           Generate all spin-off/variant data products needed...
           (full model errors, after checking resolution + manual approach)

           Try splitting 1-1.7 and 2-3 bands into smaller pieces! (low priority)
             (run `reproject_obs` and `flux_obs`, double check that I get the
              same files that Brian has sent, then make new cuts and numbers)

           For my own verification -- run full model fits w/ capped/bkg-subbed
           FWHMs, to compare to non-corrected FWHMs.  Do chi-sqr values drop
           consistently?  Does eta2 increase from ~0 values?  Possible that
           there's nothing -- but I want to check.
           (run after regions-6 is done)

* Kepler -- Get Rob/Brian to okay Kepler numbers...
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

* General note: final data for manuscript should always re-run and update
  resolution notebook

Lower priority / various checks
-------------------------------

* Manual error computation -- test it out on one/two regions, estimate how much
  time it takes (and extrapolate).  Do this after the resolution checks.

* SN 1006 -- re-run pipeline at very end and ensure all is still consistent
             re-run pipeline on _individual_ SN 1006 regions.
             Does the sub-Bohm diffusion result still stand at that level?

* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Check all constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran.  Go back to Sean's transport eq'n and rederive.

* Check transport equation for pure advection case and for damping case

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

* In plasma transport equation, could other terms matter?  E.g. the
  1/3 p df/dp dv/dx term or whatever.  How good is the simple plane approx?

Extra
-----
i.e., things that no one cares about

Write/find short scripts / hooks to git, to clear out ipynb output cells or
something... 



