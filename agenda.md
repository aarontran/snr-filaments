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

Misc questions/thoughts
-----------------------

* Can we edit our AAS abstract?
* What parameters did Sean use to compute SN 1006 profiles, or rule out damping
in SN 1006?

Main TO-DOs
-----------

In order:

2. Magnetic damping!  Run fits for all Tycho regions.
   Also fit with eta2 fixed to 1, 0.1, 10 at each scale length.
4. Use srcutlog (from Brian) to run fits with eta2 fixed at values determined
   from synchrotron break frequency.
5. Compute mE values for best loss-limited fits (to compare to data mE's).
   Do this for best damping-fits too, for completeness.
   Reason: we expect mE values for fits to be much smaller than data mE values,
   which are misleading (esp. due to large error bars on lowest energy FWHMs)

Prepare for Monday (or earlier), in order:

(NOTE: fix the automated line spectrum fit for Region 9... check what went
wrong.)

* Magnetic damping fits (chisqr as a function of `m_E`) (in progress)
* mE values from model fit curves (B-damped or loss-limited) (...)
* Downstream spectral fits w/ lines accounted for (done)
* Loss-limited fits from srcutlog (...)
* Get Brian/Rob to okay Kepler numbers (...)

(once all results are in place, go back to paper!)

Moving to lower priority:

2. Check out outliers by hand (why did the fits go wonky), and check manual
   error computation speed.
   (partial reason: B0 table is imperfect... generate new table with dense B0
   sampling and see if helps at all)
3. regenerate tables for regions-6 fits (also damping tables if need be)
4. Generate Kepler tables and run full model fits


Data processing to-dos, by individual "pipeline":
(i.e., getting FWHMs/regions/etc)

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

* Kepler -- Check full model resolution #s (ipynb)
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



