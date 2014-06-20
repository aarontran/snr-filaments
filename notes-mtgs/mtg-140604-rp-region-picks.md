Wednesday, 2014 June 4 mtg
==========================
with Rob Petre (~12:30p to 1:00p ish)

Had begun selecting regions (eyeballing) that morning, went in to take a look
and discuss.  Brought a number of questions about region picking and next
steps.

Region picking
--------------
This is for the radial profiles along Tycho's thin, nonthermal filaments.

### How to find SNR center?
Fit remnant to circle, check literature, etc...?
Didn't fully get answered, but probably just eyeball and/or fit.
See answer to following question.

### How do we get a radial distance measurement?
Motivation: figures from Ressler paper had x-axis of radial dist. What happens
if the shocks aren't aligned on a nice circle? Rob: maybe just take average,
the sensible thing to do. But, doesn't really matter -- what we need is the
distance from the shock filament; those peaks are our fiducials.  More
important to minimize curvature, minimize width.

### Does width of the region selected matter?
Not really.  Nice to be comparable, but here we really just want more signal.
(idea: if we know width, we could estimate expected curvature and hence
estimate how much widening could be due to curvature???)

### Optimizing filament region selection (interface with Python, etc)
Rob: go for it, though we definitely didn't do that for SN 1006.
Talk to Brian about this, he will probably know more.
(kind of a 2nd-order improvement)

### How far to extend region profiles, ahead/behind shock wave?
Remember, fitting only needs subsets of the profiles anyways (Fig. 7, Ressler).
Ahead of shock: we need a background subtraction, of course.
Behind shock: yes, going until you see slight increase is reasonable.

### Miscellaneous
Mentioned multiple filaments -- as expected.
Few data points, but Brian will prepare 0.5 arcsec. resolution images next week
to work with. Why can we get away without binning?  All about count rates -- in
SN1006 we had to bin for good statistics, to get a clean enough signal.
(my own note: tradeoff between SNR, and amt of signal information available?)

Random tidbit, apparently Satoru Katsuda did the region selection for SN 1006
(Ressler et al.).  Not sure if that was relevant, but just a mention.  Maybe if
we were to compare methodologies, could be useful.

Spectra and fitting
-------------------
I.e., next steps after we have these regions.

I need to get the following software:
1. specextract (bundled in CIAO), works on the .fits files to get events
2. xspec (bundled in ftools).  This requires certain calibration-like files
   from Chandra (information about effective area, spectral resolution).
   Rob says these may be generated automatically from CIAO, maybe during
   specextract stage.  Likely available from CXC website.  But anyways this is
   kind of jumping ahead.  Yes, the .fits files should have the list of
   individual events/energies/et cetera.

Spectrum fitting -- according to Rob, it does a sort of forward-fitting.
Take the model, convolve it with the detector response, obtain a chi^2 value.
Rinse and repeat. (here I need to understand the physics better)

The big step after this will be beginning to fit the profiles, for this I'll
write my own Python code or what have you.

To-dos for next few days
------------------------
(meetings, deliverables, objectives)

When region picking is done -- email Brian/Rob files, information (notes),
profiles (later this week, when you've explored this thoroughly enough).
Rob will have email access.

Play with specextract, xspec in next few days.

Email Brian about meeting Monday.  Rob will email Brian too to keep in loop.
Meet Rob/Brian both on Wednesday, get all up to speed. Rob is coming back
Monday or Tuesday, leaving Thursday next week (back to back trips to Chicago)

2nd-order effects (backburner)
------------------------------

* Optimizing region fits, again talk to Brian about this.
* (much later) after fitting profiles, try removing thermal contamination

Let's see how deep we can go (but, do it well...).
