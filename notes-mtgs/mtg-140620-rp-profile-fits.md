Meeting with Rob, morning of Friday 20 June 2014
================================================

As noted in `rsch-notes.md`, I brought in some simple profile fits.  Rob seemed
pretty happy with the look of the fits.  Some notes from discussion follow.


About fitting
-------------

## Error
Yes, it's just Poisson, sqrt(n).  Intensity scaling needs to be
accounted for, however.  I asked about the issue, that we are binning in energy
and so we have 1. different counts contribute different amts of energy, 2.
different energy uncertainties.  Yes, it probably matters -- but it's not so
important.

## Intensity units?
A.U. probably stands for arbitrary units!

## Model choice
Remember that it's strictly empirical, no physics involved.

Possibility of using interpolating spline or similar, with smoothed data?
Rob says it came up before, but they didn't use it, don't recall why.
For now, best to stick to model here to be consistent w/ previous work.

## Model choice and projection effects
On the Gaussian in SN1006 to capture projection effects -- we are not working
with spherical cows.  There is waviness, lots of structure.

Here Rob mentioned how at a conference, someone brought a window screen/mesh
kind of thing once -- look at it straight on and its transparent.
But, look at it an angle and you lose lots of light.  Now wiggle it and make
waves, and you start to see crazy patterns -- like you see in, e.g. SN1006.

So it's not unexpected that we don't need the Gaussian, here!

## Fitting routine
Rob suggests -- remember you can freeze/thaw parameters, to get your fit to
behave more nicely.


About interpretation of fits
----------------------------

## Uncertainty estimate in FWHM
Follow Ressler -- helpful to know that a delta-chi-squared of ~2.7, corresponds
to a 90% confidence interval or something similar.

## Fit quality assessment
The chi-squared will not be that great, and it's okay!  Because we just want to
pull out a width.

One issue -- the exponentials tend to shoot above the data, which may
be unphysical and would affect our FWHM calculations.

To address this, get a width out from an eyeballing approach as well -- check
that this width agrees to FWHM from fit, within error.


## Does profile height matter?
No, no, not for Ressler's model.  That's just telling you about the frequency
spectrum, in essence!  All about width.


Logistics
---------

## Satoru Katsuda
If there are questions, throw together a list for Monday, we can email Satoru!
He used to be a NASA postdoc, for some 3 years (or, was it 3 years ago?)

## Meeting
Meet with Brian and Rob on Monday -- with Nina as well.  Nina and I can learn
more about what we're each doing.  Brian needs to catch up on everything that's
happened this week, that I (we) have talked to Rob about.

