snr-research agenda
===================
Aaron Tran
Summer 2014

N.B. this is not an archival document, but rather replaces post-it notes and
other ephemeral pieces of paper.

LEARNING (2nd priority for now...)
----------------------------------

Basically an evolving to-do list.
* [snr spectra](http://www.phy.duke.edu/~kolena/snrspectra.html)
* [x-ray lines](http://www.phy.duke.edu/~kolena/strongxlines.html)

Get an intro level astrophysics textbook
1. basic radiative processes (thermal, nonthermal emission)
2. supernova physics (temperatures, radiation/spectra during evolution)
3. physics of Tycho -- interpretation of the fluff.  ISM gradient.
4. Ressler et al. (in press, ApJ)

5. X-ray telescopes.  Resolution, operation, limitations

About radiative lines (XSPEC stuff):
[atomdb](http://www.atomdb.org/Physics/units.php)

Q: why is the amplification higher than expected -- what sets the expected amt
of amplification?  (ressler mentioned expected value of 4x for unmodified,
strong shock -- where did that come from?)

(Rankine-Hugoniot conditions for strong adiabatic shock?..)

Possibly useful link on most recent galactic SNR
[here](http://chandra.harvard.edu/photo/2008/g19/media/)


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

Problem: we also have to consider spectra when optimizing

THING in TYCHO
--------------
Look at spectra of the THING...
Use a circle, for background use annulus of nearby stuff.

Disappointing.  Looks the same as the rest!  Try breaking it apart into
different pieces, see if it looks different.

About this: see Hwang and Laming (2012, ApJ 746(2)), 3rd paragraph of
introduction:

> The X-ray-emitting Si ejecta show a bipolar structure with jet-like features
> (Hwang et al. 2004; Vink 2004; Laming et al. 2006) similar to that seen in
> optical (Fesen 2001 and references therein) and infrared emission (Hines et
> al. 2004).

Agenda for code/analysis/regions
--------------------------------

4. Inspect regions + spectra/model/residuals...  Rinse and repeat

Note (June 18): I have been neglecting A LOT of parameters in XSPEC, CIAO.
Mainly in fits (e.g. fit statistic, chi vs. cstat), CIAO I haven't been
thinking about wmap, energy levels, etc... need to read more.

But, doesn't matter that much because we are just checking for no thermal
contamination, not actually using spectra in results.

Check that file hierarchy is sane.  Move whatever files need to be moved
(update scripts in time... add one or two bash scripts to pull everything into
an easy to use thing (i.e., single shot from DS9 region to output plots,
spectra, spectral fits, profile fits, figures)).


Profile fit to dos:
-------------------

7. Somewhere -- script to parse out FIT parameters / eqwidth calculations and
   make a nice formatted list to make life easier... (FOR SPECTRA)
9. get spectra for regions-4-ext
10. make regions-5 with more regions
11. Think about Tycho -- who else estimated B field amplification?
    Should the electrons be loss-limited, or what (magnetic field damping
    instead)?  Make sure our result -- decreasing filament widths -- is sane.
12. iPython notebook is kind of slow.  Refactor / clean up code


Other issue
-----------
On the matter of `m_E` computation:

My current train of thought -- keep calm and fit FWHMs.  Else, use
point-to-point geometric means (i.e., average `m_E` values, not FWHM values).

We can get `m_E` from point to point measurements (acknowledging yes there's a lot
of spread) later, but I think it's also worth doing the backtracking -- the
values of `m_E` will be more robust that way. (e.g., as seen in the
not-necessarily-meaningful global fits -- values are pretty consistent).

And, if anyone disputes the calculations of `m_E`, they can always try
reproducing it, or suggesting an alternate calculation.  Best we can do is
perform and present a simple, sensible (i.e. what-you-would-expect)
calculation.


To-dos
------

Before leaving today (monday jul 7) -- I need to have
1. manual fits to approximate model, to get some numbers (part done)
   (compare/contrast using different FWHMs for this computation)
2. rederive mE following sean's email, and email to follow up!!!! NOT DONE

Aim to send stuff to Brian tomorrow!  As I have said.
Stop hacking at the code.  Generate spectra for the regions I'm using RIGHT
NOW (forget about making new regions just yet -- maybe Sean/Steve will have
thoughts, too!)
Notate which ones can't be used with 0.7-1keV, figure out a way of doing that.


* Demonstrate that I get the same results as Sean's original `Widthfun.py`.
  This should get a cell/section to itself.

* CREATE a test case for lmfit, to verify it is doing a least squares fit as I
  would expect (just check against scipy curve fit).
  During this test -- verify that when I freeze a parameter, chi2red is
  calculated correctly with one less DOF (if it spits out a chi2red)
  This should get a cell/section.
* At some point, verify all of sean's calculated constants
  This should get a cell/section.

TUESDAY JULY 8
--------------
* See Brian's email: generate a handout of region pictures (zoomed in),
  spectra with fits and numbers, profile fit figures.  Need to send to Sean and
  Steve so they can see what we're doing.

  0.5 keV to 7 keV to verify that there is no oxygen line

* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* GENERATE SET OF REGIONS WITH GOOD 0.7-1 keV FWHM, AND SET OF REGIONS WITH BAD
  0.7-1 keV FWHM -- as discussed previously.

* In some cases, a 0.7-1 keV FWHM is found even though it's obviously a mess /
  not possible.  I need to identify and remove these regions.
  (i.e., push these to set of regions with BAD 0.7-1 keV FWHM)

* Incorporate Tycho's variable shock speed into results (Williams et al., 2013)

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
