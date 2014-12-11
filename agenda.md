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
* radio relic modeling -- lots of parallel work to look at

Shimoda, Inoue, Ohira, Yamazaki, Bamba, Vink,
submitted to ApJL: http://arxiv.org/abs/1412.2874

Main TO-DOs
-----------

* Set up extra fits for Tycho
* AAS poster
* Go through manuscript, hit the lingering TODOs and minor notes

* check expressions for cut-off energy, diffusion terms in Sean's solution
  (using Bmin vs. B0 in certain places).  Review DSA derivations, synchrotron
  derivations, assumptions, etc.

* Figure out the Bmin thing in diffusion coefficient calculations

* Eventually: review manuscript for tense, loss-limited/damping terminology /
  definitions, etc...  Check all equations, resolution, etc.
  Recompute all resolutions and update notebook

* Test the strange D^2(x) B(x) assumption? How?  Tang/Chevalier?

* The diffusion coefficient timescale depends on both upstream/downstream
  coefficients, and requires a simple integral.  What happens if we allow the
  diffusion coefficients to vary?  Does this timescale change significantly?

* Cas A -- pick subset of regions-0 for regions-1 from FWHMs/profiles/counts,
             ensuring that all regions can get decent FWHMs.
           generate new FWHM fits (1-9 keV) and get spectra of rims
           select specific regions outside rims and get spectra, need to figure
             out how much thermal contamination


### Lower priority / details
* Manuscript data (Tycho regions-6, SN 1006)
  - Tycho loss-limited fits with eta2 = 1 fixed
  - Tycho damping fits with variable distance
  - SN 1006, fits for individual regions?! (super low priority)

### Extra avenues

May not necessarily be explored.

* New EVLA data, pending reduction
* Self-similar solutions for velocity, density, B field.  Likely a marginal
  difference
* Simple hydrodynamic codes

* SN 1006 -- run pipeline on individual regions in Sean's paper

### Questions:

* Should we pitch this paper as a detailed / small scale study of Tycho's
  filaments, instead of focusing on the energy dependence result?  We try
  several approaches to constrain damping (energy-fwhm dep., spec. variation,
  radio/x-ray eyeballing), although energy-fwhm dep. is most involved.  We can
  also try exploring parameter space for spectral variation, see what
  parameters do give spectral variation consistent with observations.  But this
  is rendered less credible due to the assumption that D(x) B^2(x) is constant.

  Once this is addressed, then the intro/abstract should be rewritten
  accordingly.

* do we know anything about diffusion of GeV electrons (i.e. at radio
  wavelengths)?  What should we expect -- especially since Bohm is a lower
  limit, I'd assume it should be set by the scale/strength of magnetic
  turbulence or something...

* should we throw out the errors on our fits?  They're not very useful, now.

* how do we tie X-ray width measurements, and joint X-ray/radio profile
  modeling, together?  They are basically two separate investigations that use
  the same model/code

* Multiple filaments?  Any simple models for emission with, e.g., some simple
  sinuosoidal perturbation to a spherical shell?

* would be nice to talk about damping and relate to literature on 2-D/3-D
  simulations, more complex dynamics.  How does damping relate to plasma
  evolution, particle acceleration?

* Assuming synchrotron stuff in Tycho is predominantly from old
  shock-accelerated electrons hanging around -- does that imply the radio
  brightness should increase w/ time, as a function of the swept up mass and e-
  acceleration efficiency?  Ought to be testable, esp with the youngest
  galactic SNRs.

* could field increase?! behind the shock?! (Guo et al. 2012)

* does the strength/intensity of radio emission from SNRs require any
  amplification, or is an ambient ISM field strength enough to produce the
  observed emission?

### Radio Qs for Jack

* How badly does CLEAN mess up radio flux measurements?
* How to best use data from all configurations?  Throw all baselines together?
  (would that mess up the beam functions etc)
* NNLS algorithm to resolve shock structures?
