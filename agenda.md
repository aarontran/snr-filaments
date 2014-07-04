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
Even if we nail the undershoot/overshoot issue consistently -- we still have
the problem that our values of `$m_e$` will vary wildly from one region to the
next, a lot of it due simply to fitting uncertainty...

This gets at the issue discussed with Brian -- why the hell did they average
such disparate FWHM values, and then back out an average `$m_E$` for the
filament?

My current train of thought -- don't even worry about it!  Just fit FWHMs.  We
can get `m_E` from point to point measurements (acknowledging yes there's a lot
of spread) later, but I think it's also worth doing the backtracking -- the
values of `m_E` will be more robust that way. (e.g., as seen in the
not-necessarily-meaningful global fits -- values are pretty consistent).


Calculate individual mE values from SN 1006 to get an idea of the variability


To-dos
------
* Tabulate the results, like what Sean did in LaTeX -- put all the FWHM numbers
  together, with point-to-point `m_e` calculations and global fit `m_e`
* Send all the relevant files, plots, numbers, documentation.  Explain the
  procedure, how everything was done (for whatever paper in works).

* Profiles -- need to generate fits using intensity mosaics, because those have
  exposure and vignetting corrections.  Should make little difference though.
* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?

* Model code -- test it out first, before automating / wasting time scripting
  stuff up.  Parameter space could be weird, so might have trouble converging

* Kepler -- use CIAO `merge_obs`

Brian -- would be good to have a writeup, flow-chart like thing of what we
calculated, what we did, what equations were used where -- kind of what we're
retracing right now going through the Ressler paper.

SHOCK SPEEDS IN TYCHO ARE VARIABLE.  ACCOUNT FOR THIS IN FITS/CALCULATIONS.

* Incorporate Tycho's variable shock speed into results
* Eventually, generate a flowchart of what calculations were done, what
  equations were used, etc...
* See Brian's email: generate a handout of region pictures (zoomed in),
  spectra with fits and numbers, profile fit figures.  Need to send to Sean and
  Steve so they can see what we're doing.



