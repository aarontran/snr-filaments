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

*Work/pipeline/software engineering/whatever*

So need to
1. overhaul specextract pipeline (handle region coordinates correctly, test and
   ensure it works correctly with my adjustments)
2. rehaul data structuring for FWHMs because it's a pain in the ass to plot /
   organize
3. include switches to test different measurements/procedures, for
   reproducibility.
4. proper config files / cmd line arguments for fitting/FWHM analysis scripts?
5. clean up rsch notes document... keep track of what changes/etc I've been
   making


* overhaul pipeline for profile fit / FWHM processing again -- it's just so
  messy, and feels hard to work with

* Ah dammit.  I think I've been using inconsistent coordinate systems etc all
  this time... so my background links, spectra, regions, etc may all be
  slightly off.  Argh. Need to fix this.

* Update pipeline documentation ...

* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.
* Shift region numbering to be more logical (instead of 1, 10-13).

* Eventually: run whole pipeline on one set of ALL regions, sampled all around
  SNR, save the output.  Use this to argue/show why regions are good/bad.

Model stuff
-----------

One major issue -- appropriate electron distribution parameters for TYCHO!
Cut-off energy, spectrum index, etc.  Where did he get value 2.2 for SN 1006?

One big question -- what are the qualitative differences between the models --
this seems somewhat addressed by the discussion of electron cut-off energy.

Why should the cut-off energy be set by equating loss times and acceleration
times?  Is that to say that, the electrons radiate faster than they can be
accelerated?  But that seems weird -- because I thought the acceleration is
what gives rise to the radiation.  Please read rettig/pohl carefully, take
notes, rederive things.

How the hell to get such strong B field values?!

Another issue -- when converting from particle energy E to radiated frequency
nu, doesn't the use of `$\nu_m = c_m E^2 B$` implicitly invoke the
delta-function approximation for synchrotron radiation?  I have no idea how
you'd do it otherwise though -- the diffusion coefficient is set by a singular
electron energy, which corresponds to some distribution of synchrotron
radiation.  But, perhaps, for an individual electron this is not a bad
approximation?  I don't understand the physics of the derivation of synchrotron
stuff (MUST REVIEW ALL OF THIS...).

Further clarification with Sean -- after having rederived `m_E` and all, I am
convinced that the output from `Widthfun.py` should allow us to use equation
(23) to compute `m_E`.  I think, Sean misread my question / my question was
unclear, and so that is important in attempting to describe `m_E` using the
FULL numerical code (the full power of this orbiting ........)

### Python lmfit, approx equation

* CREATE a test case for lmfit, to verify it is doing a least squares fit as I
  would expect (just check against scipy curve fit).
  During this test -- verify that when I freeze a parameter, chi2red is
  calculated correctly with one less DOF (if it spits out a chi2red)
* At some point, verify all of sean's calculated constants
* Demonstrate that I get the same results as Sean's original `Widthfun.py`.

### Full model

* Check transport equation for pure advection case
* Check numerical prefactor 8.3 TeV for electron cutoff energy
* Update numbers from Pacholczyk (?), consider adding more entries (can ask
  Sean about this)
* Check FORTRAN code for memory-efficient array indexing?
* Consider caching / memoizing tabulated electron distributions?
  Problem: advective/diffusive lengthscales depend on B0, mu, eta, etc.
  So, as we try new parameters, the electron distributions change...
  which does physically make sense, after all.  sigh.
  How can we get around the speed issue?

### All models
* Incorporate Tycho's variable shock speed into results (Williams et al., 2013)
* USE 3 kpc instead of 2.3 kpc, and will need to investigate this effect soon
  (2.3 to 4 kpc range).  Brian's paper actually favored a larger distance.

* Model code -- test it out first, before automating / wasting time scripting
  stuff up.  Parameter space could be weird, so might have trouble converging

Then look at spatial dependence of B field, robustness of numbers from approx
python model. Brian asked, does B field scale with stronger energy dependence
Should verify this. Need to have discussion/commentary on what affects what.

On using different calculations of the FWHM: how do I show the effect of these
different procedures?  Some kind of normalization? (Figure 10 of Ressler).
Quantify effects on calculation of `m_E`, B0, eta2.


Documentation/poster/paper
--------------------------

Start making poster later this week, have nice draft ready by Monday
(print by Wednesday morning at latest -- setup 2-6p Weds, 7:30-11a Thurs)
(presentations: 11a-2p, must be present 12a-2p)
(spec: 45" by 45" max)

More supernova remnants
-----------------------

* Kepler -- use CIAO `merge_obs`
* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me... (so, verify that it works right with Tycho, then
  test on kepler???)
* Question: looks like `reproject_obs` will give evt file (which can then be
  partitioned by `dmcopy`, is that good enough....

Read `merge_obs` documentation, play with outputs
