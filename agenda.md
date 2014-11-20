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

Speculation:
* effect of varying electron population spectrum (besides Zirakashvili and
  Aharonian, see also new arXiv paper by Pohl, Wilhelm, Telezhinsky).
  http://arxiv.org/abs/1411.2891 (figures are missing...)
* radio relic modeling -- lots of parallel-seeming work to look at!...

Main TO-DOs
-----------

Keep end game in mind, where do we cut it off and start a new paper?

* WRITEUP results, edit manuscript w/ fixes etc
* Check spectral variation from models (damping case for R&P)

* Radio reduction...
* Test the strange D^2(x) B(x) assumption
* Simple hydrodynamic codes?
* Self-similar solutions for veloc, B behind the shock,
  to get radio/X-ray profiles farther downstream

* check expressions for cut-off energy, diffusion terms in Sean's solution
  (using Bmin vs. B0 in certain places).

### Lower priority / details
* Manuscript data (Tycho regions-6, SN 1006)
  - Tycho loss-limited fits with eta2 = 1 fixed
  - Tycho damping fits with variable distance
  - SN 1006, fits for individual regions?! (super low priority)
* Add scale bar (arcsec) to Tycho SNR image...
* Check outliers by hand (why did the fits go wonky)
* ODE solving, full model with damped B-field (see Tang/Chevalier)

### Radio Qs for Jack

* How badly does CLEAN mess up radio flux measurements?
* How to best use data from all configurations?  Throw all baselines together?
  (would that mess up the beam functions etc)

### Theory stuff

* The diffusion coefficient timescale depends on both upstream/downstream
  coefficients, and requires a simple integral.  What happens if we allow the
  diffusion coefficients to vary?  Does this timescale change significantly?
* Rederive the advective (transport) solution in damped and non-damped cases...

* Solution to transport equations -- Tang and Chevalier (arxiv:1410.7510v1,
  posted October 28) compute a much more general solution (time-dependent,
  arbitrary injection spectrum, using an inverse Laplace transform computed
  numerically with Talbot's method?!!! (geez)

* Diffusion in radio?

* General note: final data for manuscript should always re-run and update
  resolution notebook

Lower priority / various checks
-------------------------------

* Cas A -- pick subset of regions-0 for regions-1 from FWHMs/profiles/counts,
             ensuring that all regions can get decent FWHMs.
           generate new FWHM fits (1-9 keV) and get spectra of rims
           select specific regions outside rims and get spectra, need to figure
             out how much thermal contamination

* SN 1006 -- re-run pipeline at very end and ensure all is still consistent
             re-run pipeline on _individual_ SN 1006 regions.
             Does the sub-Bohm diffusion result still stand at that level?

* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Rederive synchrotron spectrum results from scratch (carefully accounting for
  factors of sin theta).  All of Sean's synchrotron constants are confirmed
  (although not better than 0.1% as Pacholczyk's values are old + due to
  roundoff, etc).

* Check transport equation for pure advection case and for damping case

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

For data storage / archiving / keeping track of variant computations -- learn
to use SQL databases.  Or pandas dataframes.

