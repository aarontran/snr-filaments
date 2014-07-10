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

(Rankine-Hugoniot conditions for strong adiabatic shock?..)

turbulent B field amplification (nonlinear plasma physics):
doi:10.1038/nphys2978 (in a lab! this is really cool)
[ppt](http://fermi.gsfc.nasa.gov/science/mtgs/symposia/2007/p4/P4.1_Ellison.pdf)
[AR Bell](http://mnras.oxfordjournals.org/content/353/2/550.full.pdf) (highly
cited paper on this stuff)

Possibly useful link on most recent galactic SNR
[here](http://chandra.harvard.edu/photo/2008/g19/media/)

Organize literature collection/review eventually...


Optimization of region selection
--------------------------------

2nd order stuff (nice but less important)
* Optimization
* Calculate *expected* curvature based on width of regions
  (idea, could quantify how wide these may be. But, multiple filaments likely
  to be a larger confounding effect)

Idea: try testing projections for minimum peak width, e.g. using pyds9.
Could find when peak height is maximized, when trough behind peak is minimized
(ratio of trough/peak height), estimate peak width, etc...

Supply initial guess regions, then rotate them about their centers (?) until
the filament peak width is minimized.  You should generate the projections in a
consistent way (just for ease of use/manipulation), e.g. starting from outside
and moving in.

Problem: we have to consider spectra + results in 5 bands

THING in TYCHO
--------------

Try breaking it apart into different pieces, see if it looks different.

About this: see Hwang and Laming (2012, ApJ 746(2)), 3rd paragraph of
introduction:

> The X-ray-emitting Si ejecta show a bipolar structure with jet-like features
> (Hwang et al. 2004; Vink 2004; Laming et al. 2006) similar to that seen in
> optical (Fesen 2001 and references therein) and infrared emission (Hines et
> al. 2004).




Stuff to send around
--------------------

* NEED IMAGES OF TYCHO WITH REGIONS (zoomed in and labeled)

* we can put the plots on a website instead, doesn't have to be a PDF doc.
  Something for collaborators only, or what have you.  nbviewer looks okay.

* Yes, go ahead and make your region changes / whatnot before sending.
* Make tickmarks readable, fiddle with formatting/sizing
* Calculate spectra with all 750 ks of data.  Ask Nina for her modified copy
  of Brian's script, which she got to work on a mac.
* FWHM comparison figure?
* Add tables of FWHM values, for people to see.  Eventually it will have to be
  in LaTeX.

* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.
* Shift region numbering to be more logical (instead of 1, 10-13).

At some point, walk through EACH region with ALL people
and verify that selections are good.


Model stuff
-----------
2. rederive mE following sean's email, and email to follow up!!!! NOT DONE

* CREATE a test case for lmfit, to verify it is doing a least squares fit as I
  would expect (just check against scipy curve fit).
  During this test -- verify that when I freeze a parameter, chi2red is
  calculated correctly with one less DOF (if it spits out a chi2red)
  This should get a cell/section.
* At some point, verify all of sean's calculated constants
  This should get a cell/section.

* Demonstrate that I get the same results as Sean's original `Widthfun.py`.
  This should get a cell/section to itself.

* Incorporate Tycho's variable shock speed into results (Williams et al., 2013)

* USE 3 kpc instead of 2.3 kpc, and will need to investigate this effect soon
  (2.3 to 4 kpc range).  Brian's paper actually favored a larger distance.

* QUESTION: on computing `m_E` point-to-point, why use lower energies, instead
  of say using median energy of band (median energy of photon population -- of
  course that will be pretty close to lower energy, I guess, esp. in the
  farther spectrum)

AFTERWARDS
----------

* Tabulate the results, like what Sean did in LaTeX -- put all the FWHM numbers
  together, with point-to-point `m_e` calculations and global fit `m_e`.
  Send plots, numbers, documentation -- explain procedure, how everything was
  done.
* Eventually, generate a flowchart of what calculations were done, what
  equations were used (where used), etc... kind of what we're retracing now
  going through the Ressler paper.

* Model code -- test it out first, before automating / wasting time scripting
  stuff up.  Parameter space could be weird, so might have trouble converging

* Kepler -- use CIAO `merge_obs`
* sanity check that when I run CIAO `merge_obs`, I get the same files that
  Brian has been sending me...

Later -- integrate region names/labels into workflow?

On using different calculations of the FWHM: how do I show the effect of these
different procedures?  Some kind of normalization? (Figure 10 of Ressler).
Quantify effects on calculation of `m_E`, B0, eta2.

After sending spectra to brian, then look at spatial dependence of B field,
robustness of numbers.  Brian asked -- does B field scale with stronger energy
dependence?  Should verify this.
Need to have discussion/commentary on what affects what.
