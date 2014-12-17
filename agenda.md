snr-research agenda
===================
Aaron Tran

High level questions to consider
--------------------------------

From poster session (July 31):
* Ori -- sanity check on checking luminosity and shock kinetic energy, though
  Rob/Brian note that the synchrotron radiation is very inefficient (esp.
  compared to when the remnant cools, H starts recombining and emitting thermal
  lines, T ~ 10^4 K)

Other stuff:
* effect of varying electron population spectrum (besides Zirakashvili and
  Aharonian, see also new arXiv paper by Pohl, Wilhelm, Telezhinsky).
  http://arxiv.org/abs/1411.2891 (figures are missing...)
* radio relic modeling -- lots of parallel work to look at
* CR accel efficiency: Shimoda et al., submitted to ApJL, arxiv:1412.2874

Main TO-DOs
-----------

### Regular to-dos

* check expressions for cut-off energy, diffusion terms in model solutions
  (using Bmin vs. B0 in certain places).

* Keep reviewing manuscript, TODOs/prose

* The diffusion coefficient timescale depends on both upstream/downstream
  coefficients, and requires a simple integral.  What happens if we allow the
  diffusion coefficients to vary in space?  Does this timescale change significantly?
  (this requires a DSA re-derivation)

* Test the strange D^2(x) B(x) assumption? How?  Tang/Chevalier?

* Eventually: review manuscript for tense, loss-limited/damping terminology /
  definitions, etc...  Check all equations, resolution, etc.
  Recompute all resolutions and update notebook

* Cas A -- pick subset of regions-0 for regions-1 from FWHMs/profiles/counts,
  ensuring that all regions can get decent FWHMs.
  Generate new FWHM fits (1-9 keV) and get spectra of rims
  Select specific regions outside rims and get spectra to estimate amt of
  thermal contamination (reflected light etc)
  At this point, just for fun.  May not go anywhere.

### Telecon things

* How to address change to SN 1006 results?  Currently 1 paragraph discussing
  how width-energy doesn't discriminate here.

* __poster/paper organization__: could we pitch paper as a detailed / small
  scale study of Tycho's filaments, instead of focusing on energy dependence?
  We try several approaches to constrain damping (energy-fwhm dep., spec.
  variation, radio/x-ray eyeballing), although energy-fwhm dep. is most
  involved.  Once addressed (or not addressed), rewrite intro/abstract.

* do we know anything about diffusion of GeV electrons (i.e. at radio
  wavelengths)?  What should we expect -- especially since Bohm is a lower
  limit, I'd assume it should be set by the scale/strength of magnetic
  turbulence or something... no idea honestly.

* should we throw out the errors on our fits?  They're not very useful, now.

* Check w/ steve about VLA A config (any merging of multiple configs / freqs?)

* Ideas for (distant) future studies or work?

* Testing the strange D^2(x) B(x) assumption?

### General Qs/things

Run by Brian...

* Intensity comparison from modeling.  Haven't done this.

* Hard drive for mac to copy all my files / material off?

* Run plots of FWHM fit profiles + data by Brian (+Rob) briefly.

* Kepler radio images?  Something I can do very quickly.
  2004/2005 data (PI: Delaney?)?  is it published?
  http://adsabs.harvard.edu/abs/2009AAS...21348804G ...

* Assuming synchrotron stuff in Tycho is predominantly from old
  shock-accelerated electrons hanging around -- does that imply the radio
  brightness should increase w/ time, as a function of the swept up mass and e-
  acceleration efficiency?  Ought to be testable, esp with the youngest
  galactic SNRs.

Major question I have -- what would we expect for shock evolution, dynamical
age of shock, etc?  I am wondering if there is a correlation between:
- loss-limited rim (radio rise)
- stronger width-energy dependence, perhaps comparatively dimmer in X-ray vs.
damped rims
- dynamically young shock, or shock that hasn't interacted with ejecta /
other material?
look at the really young / faint shocks... it's real hard to say for the
fainter ones.  There's sometimes simply no rim in soft X-ray to be seen.

This is basically based on two observations:
1. Regions J, L, M around WNW of Tycho.  J is associated w/ regions 11,12
   in X-ray, which have harder spectra than rest of Tycho (smaller gamma,
   larger cutoff energy).  Both regions (12 especially) require the
   "weak-field" damping in order to get correct energy dep.  Fit is _not
   physically meaningful_ (weak-field damping not consistent w/ radio rise),
   but suggests stronger energy dependence than in other regions.
   Region J is pretty bright, next to ejecta knot, though.  L, M are fainter
   and don't show clear soft X-ray rims.
2. Region O (or, X-ray Region 18).  This is the faint, fast moving shock (maybe
   1.5x to 2x faster than similar shock to southeast).  Again, strong energy
   dependence in X-ray (weak-field damping fit), here there was a (faint) rim
   though.  Definitely weaker than filament rims rather SE.

All extremely speculative but I would be curious to know how damping plays into
a shock's evolution.  Do turbulence cascades occur in simulations?

### Broader exploration

Need to learn more (basic undergrad/early grad astrophys.) first.

* General question for myself: what can optical observation tell us about
  shocks?  Balmer lines -> shock velocity (broadening?)

* How does ejecta alter nonthermal radiation -- field lines?  Here see Guo+
  (2012).  They find RT instabilities at CD give very strong magnetic fields
  (this, of course, makes sense -- very bright regions associated w/ ejecta at
  FS).  Ties into my question about shock ages...

* Does the strength/intensity of radio emission from SNRs require any
  amplification, or is the advected, initially shock-compressed B field enough
  to produce the observed emission in the remnant interior?

* relate damping to literature on 2-D/3-D simulations, more complex modeling.
  How does damping relate to plasma evolution, particle acceleration?

* Multiple filaments?  Any simple models for emission with, e.g., some simple
  sinuosoidal perturbation to a spherical shell?

### Extra avenues

May not necessarily be explored.

* New EVLA data, pending reduction
* Self-similar solutions for velocity, density, B field.  Likely a marginal
  difference
* Simple hydrodynamic codes
* SN 1006 -- run pipeline on individual regions in Sean's paper

### Radio Qs for Jack

* How badly does CLEAN mess up radio flux measurements?
* How to best use data from all configurations?  Throw all baselines together?
  (would that mess up the beam functions etc)
* NNLS algorithm to resolve shock structures?
