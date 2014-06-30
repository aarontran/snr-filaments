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

Get an intro level astrophysics textbook `w__w`.
1. basic radiative processes (thermal, nonthermal emission)
2. supernova physics (temperatures, radiation/spectra during evolution)
3. physics of Tycho -- interpretation of the fluff.  ISM gradient.
4. Ressler et al. (in press, ApJ)

5. X-ray telescopes.  Resolution, operation, limitations

About radiative lines (XSPEC stuff):
[atomdb](http://www.atomdb.org/Physics/units.php)

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

(remember that we get THREE projection plots out, one for each band...)

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
thinking about wmap, energy levels, etc... need to read more `@_@`.

But, doesn't matter that much because we are just checking for no thermal
contamination, not actually using spectra in results.

Check that file hierarchy is sane.  Move whatever files need to be moved
(update scripts in time... add one or two bash scripts to pull everything into
an easy to use thing (i.e., single shot from DS9 region to output plots,
spectra, spectral fits, profile fits, figures)).

Profile fit to dos:
-------------------

3. Verify that other functions give results within error bars
   (at least, spline / eyeball fit of sorts)
   (this is a major concern -- many cases of peaks overshooting or
   undershooting the actual data)

4. Generate count images for 1-1.7 keV instead of 1-2 keV.
   Look at the profiles and estimate how many counts are lost.
   If that's workable, we can just use that band and ignore the Si issues.
   (can you quantify this, give some numbers?)

6. Twiddle with fits / robustness / etc more??!???

7. Somewhere -- script to parse out FIT parameters / eqwidth calculations and
   make a nice formatted list to make life easier... (FOR SPECTRA)

Idea
----
What if we try the binned data -- will that help the uncertainties???

Spline fourier filtering: http://bigwww.epfl.ch/publications/unser9301.pdf

Issue
-----
Really, really need to address this issue of overshooting/undershooting the
peak, because it will add more uncertainty, than the stretching procedure will.

Now that the uncertainties are suddenly nailed down -- this brings us to a host
of other problems, namely the 1. questionable fitting procedure, 2. ignorance
of the effect of varying peak height (stretch up/down)

Some possible approaches
1. make your damn fitting routine better
2. do the stretch thing, but in the y-axis now.
