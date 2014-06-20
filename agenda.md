snr-research agenda
===================

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

THING in TYCHO
--------------
Look at spectra of the THING...
Use a circle, for background use annulus of nearby stuff.

Disappointing.  Looks the same as the rest!  Try breaking it apart into
different pieces, see if it looks different.

Agenda for code/analysis/regions
--------------------------------

Version control all this stuff (wait until transferred to imac)

4. Inspect regions + spectra/model/residuals...  Rinse and repeat

Note (June 18): I have been neglecting A LOT of parameters in XSPEC, CIAO.
Mainly in fits (e.g. fit statistic, chi vs. cstat), CIAO I haven't been
thinking about wmap, energy levels, etc... need to read more `@_@`.


MADE NEW BACKGROUND SPECTRA -- BUT, NEED TO RELINK SPECTRA AND CHECK THAT
THINGS LOOK OKAY.

Check that file hierarchy is sane.  Move whatever files need to be moved
(update scripts in time... make notes here in the agenda).

Profile fit to dos:
-------------------
2. Twiddle with fit, test robustness, try other parameters/models?
3. Calculate FWHMs and errors on FWHMs (how does this stretch thing work)
4. Consider other approaches to estimate FWHM.  Try eyeball approach, try
   spline.  Ensure calculated FWHMs fall within est. error of our model FWHM.
   (Rob mentioned spline/interpolation/similar had been mentioned before
   not sure why not used? but fit to be consistent with prev. paper)

Compile questions for Satoru Katsuda (former postdoc..)

run regions-good-4.  Twiddle the regions with an eye towards getting more
counts and getting cleaner profiles (while maintaining spectrum cleanliness)
