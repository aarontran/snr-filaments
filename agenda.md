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
(Drury mentions compression ratio @ shock in DSA...)

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


Thought -- should we be constraining the resolution... checking that it
works/makes sense?  How do we quantify when the grid is sufficiently resolved?!
Compute the error from changing the grid fractionally (e.g. doubling)... kind
of like in Math 228B...

To be clear -- Sean used ecut in all of the model fits?

How do you rule out magnetic damping, exactly?  Looking at Figure 4, with ab =
0.05 rs, it seems like that would give you a decent size mE.  How did Sean get
the numbers that he did (i.e. that mE must be order -0.1 for damping)?

Question: is mu restricted to fall within `[0, 2]` (as is suggested by the
turbulence - mu relation given in Sean's paper?).

In table 8 -- how do you define unobtainable?
e.g., for mu =0.5 I can get a chisquared of at best 6.2 (compared to abt 4 for
the higher values of mu). (with manual fitting, that is)

One major issue -- appropriate electron distribution parameters for TYCHO!
Cut-off energy, spectrum index, etc.  Where did he get value 2.2 for SN 1006?
Answered by Sean's email -- see NRAO ERA (power law e- begets power law
synchrotron... work backwards. but there's also zirakashvili and aharonian
(2007) to get more detail?)

Another question: how do you know the compression ratio is 4?  It seems like a
nice assumption -- what happens if we change that assumption? where does that
factor into the model? (larger compression ratio -- smaller plasma velocity --
shorter advective lengthscale -- could help explain narrow filaments, although
not the scaling!)  (see brian's paper...)

Solution: throw this into the simple model... some numbers:
> compression ratio = 4
> Filament 1: mu = 1.00    eta2 = 1.06 +/- 0.37    B0 = 112.55 +/- 4.27
> Filament 2: mu = 1.00    eta2 = 0.08 +/- 1.05    B0 = 144.84 +/- 26.26
> Filament 3: mu = 1.00    eta2 = 0.01 +/- 0.14    B0 = 80.88 +/- 2.58
> Filament 4: mu = 1.00    eta2 = 0.01 +/- 0.06    B0 = 118.17 +/- 2.40
> Filament 5: mu = 1.00    eta2 = 4.81 +/- 5.18    B0 = 155.77 +/- 36.22
> 
> compression ratio = 6
> Filament 1: mu = 1.00    eta2 = 0.47 +/- 0.16    B0 = 85.89 +/- 3.26
> Filament 2: mu = 1.00    eta2 = 0.03 +/- 0.47    B0 = 110.53 +/- 20.05
> Filament 3: mu = 1.00    eta2 = 0.01 +/- 0.06    B0 = 61.89 +/- 2.07
> Filament 4: mu = 1.00    eta2 = 0.01 +/- 0.00    B0 = 90.42 +/- 2.22
> Filament 5: mu = 1.00    eta2 = 2.14 +/- 2.30    B0 = 118.88 +/- 27.64
> 
> compression ratio = 8
> Filament 1: mu = 1.00    eta2 = 0.27 +/- 0.09    B0 = 70.90 +/- 2.69
> Filament 2: mu = 1.00    eta2 = 0.02 +/- 0.26    B0 = 91.25 +/- 16.54
> Filament 3: mu = 1.00    eta2 = 0.01 +/- 0.03    B0 = 51.28 +/- 1.82
> Filament 4: mu = 1.00    eta2 = 0.01 +/- 0.06    B0 = 74.92 +/- 2.12
> Filament 5: mu = 1.00    eta2 = 1.20 +/- 1.30    B0 = 98.13 +/- 22.82
> 
> compression ratio = 20
> Filament 1: mu = 1.00    eta2 = 0.04 +/- 0.01    B0 = 38.49 +/- 1.46
> Filament 2: mu = 1.00    eta2 = 0.01 +/- 0.01    B0 = 50.91 +/- 10.77
> Filament 3: mu = 1.00    eta2 = 0.01 +/- 0.00    B0 = 29.00 +/- 2.27
> Filament 4: mu = 1.00    eta2 = 0.01 +/- 0.00    B0 = 42.36 +/- 3.41
> Filament 5: mu = 1.00    eta2 = 0.19 +/- 0.21    B0 = 53.27 +/- 12.39

Similar -- empirically check effect of changing remnant distance (2 to 4 kpc),
how do FWHMs and other parameters, change w/ distance?  Our guess/expectation:
scaling shouldn't change, just magnetic field estimate.  As Brian said use 3
kpc for now, scale his shock speeds appropriately.

Check all your constants urgh.  `snr_catalog.py`, model fitting code both
wrapper and fortran

Update code deep review eventually (correction to the negative sign in electron
distribution functions, explain the Ecut matter, etc)

One big question -- what are the qualitative differences between the models --
one aspect is addressed by the discussion of electron cut-off energy.

Why should the cut-off energy be set by equating loss times and acceleration
times?  Is that to say that, the electrons radiate faster than they can be
accelerated?  But that seems weird -- because I thought the acceleration is
what gives rise to the radiation.  Please read parizot et al (2006), take
notes, rederive things.

How the hell to get such strong B field values?!
Solution: read the damn literature

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
