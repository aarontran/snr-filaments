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

Learning (practical stuff)
--------------------------

Maybe useful...
[GSL manual, section on nonlinear fitting](https://www.gnu.org/software/gsl/
manual/html_node/Nonlinear-Least_002dSquares-Fitting.html)


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

### Software engineering blargh

1. overhaul specextract pipeline (handle region coordinates correctly, test and
   ensure it works correctly with my adjustments)
2. rehaul data structuring for FWHMs because it's a pain to plot / organize
3. include switches to test different measurements/procedures, for
   reproducibility.
4. proper config files / cmd line arguments for fitting/FWHM analysis scripts?
5. clean up rsch notes document... keep track of what changes/etc I've been
   making

* overhaul pipeline for profile fit / FWHM processing again -- it's just so
  messy, and feels hard to work with

* I think I've been using inconsistent coordinate systems etc... background
  links, spectra, regions, etc may all be slightly off (physical coords may
  differ w/ single obsID vs. merged obsID...).  Can't remember the exact issue,
  need to double check.

* Modularize code (separate out components for profile fitting, spectra
  extraction, region processing, whatever).  Make it easier to keep track of.

* Update pipeline documentation ...

### Actual science

* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.
* Shift region numbering to be more logical (instead of 1, 10-13).
* Slice smaller regions

* Eventually: run whole pipeline on one set of ALL regions, sampled all around
  SNR, save the output.  Use this to argue/show why regions are good/bad.

On using different calculations of the FWHM: how do I show the effect of these
different procedures?  Some kind of normalization? (Figure 10 of Ressler).
Quantify effects on calculation of `m_E`, B0, eta2.

Bibdesk -- mark up all the papers I've really gone through and understood.
Then, add the relevant citations and notes, to be sure I haven't missed or
misunderstood anything.


Models for filament widths
--------------------------

### Main / high-level TO-DOs

* Start geometric average FWHM computation, w/ appropriate errors
* Start running Kepler pipeline (will help organize code, refresh mind)
  Data structure/storage / pipeline review, how to store output FWHM data?
  Clean up current data in limbo
* When tired/sleepy, make big tables for appendix... and fix rounding/precision

* send new copy to Rob/Brian whenever/if ready, before monday...
* Bring some kepler data/images on monday...


### Weekend plan

* Keep generating tables, and putting together tables into paper
  (some tonight)
* Make plots of geometric FWHMs, errors, and fits
* Finish writing through rest of methods and results, send to Rob/Brian
  BEFORE SUNDAY.

* Review pipeline layout, code, READMEs.  Draw it out on paper.
  **address coordinate system inconsistency in ds9, check mergeobs docs**
  Write the code to save fitting data.  Write code to generate and merge tables
  automatically (FWHM tables from Fitter object; best fit tables of
  parameters/errors).  Human intervention only for final cleanup etc.

* Start running Kepler pipeline Saturday night.  Specextract needs several
  hours, and there will probably be minor issues associated w/ adapting code to
  Kepler


### List of plots/subplots/tables
See my email to Rob/Brian!

Main text:
* table, geometric average FWHMS + errs + m\_E for each filament
  *done, needs to be formatted*
* table, best fits for each filament
  *computing, 2014 Sep 13*
* plots, best fits for each filament
  *computing, 2014 Sep 13*
* (table, average over flmts of best fits for each region)
  *need flmts 2-5*
* (table, best fit for each filament with eta2 = 1)
  *done, LaTeXed*
* (table, best simple model fits for each filament -- highlight differences)
* (table, best fit for global average)
* (plot, best fits for global average of FWHMs)

Supplement/appendix:
* table, profile widths & chi-squared & m\_E for each region
  *done but needs to be formatted*
* plots, profile widths/fits for each region
  *rerun regions 11-13*
* table, spectral fit params + line detection/width for each region
* (table, best fit parameters for each region)
  *rerun regions 11-13*
* (table, plots, best fits w/ different averaging method)
  *need FWHM table; fits, plots done (arithmetic averaging)*
* (table, best fit parameters w/ different FWHM calculation method)
* (SN 1006 comparison to Sean's data)
  *rerun SN1006, flmt 4, 3 energy bands for completeness*


### Paper/table to-dos
* LATEX TABLE GENERATING -- List number of DOFs where relevant.
  or report chisqr red directly.
* NOTE my automated rounding is kinda sucking... ugh.  Trust the first 2 or 3
  digits (3 for numbers, 2 for errors), everything else is suspect.
  Fix this...

* maybe useful to compute advection/diffusion lengths from fits?
  Better yet -- ratio of advection/diffusion lengthscales
  problem is that this IS energy dependent.

* ALSO we need `m_E` from this, and the data...
  (Rob, Brian: report this in the results, and make it clear
   what's going on...
   I also want to report `m_E` as estimated from MODEL fits, which to be
   clear is a separate "observable" that I expect to be more robust
   than the point-to-point measurements)

  Also, print out `m_E` values, point to point and from `width_dump` model...
  see some of the old notebooks in `code-profiles/` (merge functionality
  together and throw out old / unused things)


### Code to do this better

I need a data structure to store best fit parameters + errors if available.
Critical for data access, manipulation, formatting (stitching tables together
nicely), checking precision, further python operations, etc.

In particular, separate data structure allows us to split data generation and
plot/table generation completely.  Keep the methods to spew out stuff in
iPython notebook, for interactive usage / immediate checking while generating.
(at least, for now, until we can get plotting/stuff working nicely -- once
the data is saved to a usable format then it's easy to separate)

How to best store intermediate numbers?  Some options I'm considering:
1. copy paste into plain text file
   need to include some critical metadata:
   * what pre-generated table was used?
   * what SNR settings (numbers etc) were used?
   * what error keywords / args were used?
   * which commit should we go back to?
   * have I modified the data from its original/raw state?
2. store pkl of numpy array (The one that gets fed into tables), with
   string table output, etc.  Also needs metadata (as above)
   But, easier to manipulate, plot again, etc.
3. spreadsheet of data (+ metadata of course)

Later, more stuff for pipeline:
* routine to manually verify data and errors...
* routine to save config files (for fitting numbers...)
* routine to merge tables together


### VARIOUS CORRECTIONS, QUESTIONS, ETC

* FWHM errors must be 68.3% CI errors, not 90% CI, for covariance matrix /
  chisqr from fits to individual regions.
  ALSO, for nonlinear fit we are not guaranteed a correspondence between
  Delta-chi-squared and CI errors; I note that Satoru didn't mention a
  confidence interval % explicitly.  So better just to say
  Delta-chi-squared = 1.
* Averaging FWHMS -- must take error following Student's t, not just stderr
  Same problem, actually, if averaging fitted `B0`, `eta_2` values.
  Need to fix that too...
* Does it matter that we are placing the data points at the BOTTOM of the
  energy band?  Intuitively makes sense, esp. as energies fall off.
  I remember discussing this w/Brian at some point.
  But, would the 0.7-1 keV and 1-1.7 keV bands be better centered elsewhere?
  Especially relevant for fitting.  But maybe better to just be consistent.
* Thought -- in regards to the overshooting matter -- yes, it could overshoot,
  but the total integrated flux should be the same.  I think that tells us,
  there is a limit to how much it can overshoot. (FOR PROFILE FITTING)

Key question underlying all of my error fretting:
what's the correct way to think about the statistics and measurements, so that
we can interpret our model fit results correctly?...


### Various checks

* Full model code -- check resolution error due to Pacholczyk tables (!)
  Might also be worth double checking some of the internal numerical integrals.
  And, check how accurate/correct Pacholczyk's Bessel function values are now.
  (ASK SEAN -- ANY REASON FOR TABLE SELECTIONS?)
* Ask Keith Arnaud about error statistics... but, lower priority
  At what chi-sqr does our simple confidence interval analysis break down/fail?
  I don't understand the theory behind this very well.
* Remnant radio spectral indices?!
  See [ppt slides](http://www.astro.le.ac.uk/~cbp1/cta/Talks/TonyBell.pdf),
  just pick some numbers and report in text (saave in SNR catalog, w/ sources
  and explanation for choices).
* Check all constants.  `snr_catalog.py`, model fitting code both
  wrapper and fortran.  Go back to Sean's transport eq'n and rederive.
* Update code deep review eventually (discuss: correction to the negative sign
  in electron distribution functions, explain the Ecut scaling / calculation)

* Check transport equation for pure advection case
* Check numerical prefactor 8.3 TeV for electron cutoff energy

* Check: sean used energy cutoff in all his model fits? (table 8)


### General (higher level to-dos)

* Remember brian's suggestion (from Friday july 25): how does mE depend on
  energy? what happens if you fit a straight power law to that???

* Look at azimuthal dependence of B field, robustness of numbers from approx
  python model. Brian asked, does B field scale with stronger energy
  dependence? (NOT ADDRESSED as of July 30)
  (answer: yes, basically, for roughly constant FWHMs. more energy dependence
  requires more B field, less eta2)


### Kepler, CIAO stuff

* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me...
* Question: looks like `reproject_obs` will give evt file (which can then be
  partitioned by `dmcopy`, is that good enough....
* Read `merge_obs` documentation, play with outputs


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
* How to get larger compression ratios?  What does it imply?
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



