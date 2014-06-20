Wednesday, 2014 June 18 mtg
---------------------------
with Rob Petre (~11:50a to 12:30p ish)

Brought in profile plots, spectra for `good-3` regions to review.
A few salient remarks:

Region spectra, thermal emission
--------------------------------

* The fact that the southern filaments are mushed out is interesting!  Esp. if
  it is just synchrotron emission.  Tells us something is going on...
* We will need to show that we tried these regions, found nothing!  Explain
  what we found...

I noted the Si/S line contamination in the spectra, although that is very weak.
Rob says we will have to estimate/bound the thermal flux... a few ways to
attack this issue:

1. plot spatial profiles for silicon photons.  Can try to show that, e.g., it
   is confined to the back region and doesn't affect our synchrotron rims.
   So we'd expect the flux of Si photons to peak further back, away from the
   main energy rims.
2. excise the counts from the line regions -- ask Sean to modify his model
   appropriately, to account for missing counts... (then we'd use different
   energy bands, essentially?).  Compare to results with the small lines.
3. general approach to estimating error is as follows:
   first we have the basic fit.
   try fitting with thermal component, see if that improves the fit
   then, add artificial thermal components until the multi-component fit is
   superior / is clearly a better model match... (look at chi2)
   then we can say, okay, here's the maximum thermal flux that could be present
   in our data, but we know it's lower because our fit quality only improves so
   much...

Rob says, he thinks it will probably be a wash.  Whatever way we do it, we'll
get some kind of small error, which (if we are doing it right) should be
comparable across all these approaches.

Since as Rob said, Tycho is ~50% silicon photons... if we can bound the silicon
emission to be small, the rest is probably small too -- very good for us!


Profile fitting and more stuff
------------------------------

Rob: start with the really good profiles for fitting, will spend some time on
this.

Q: How do we get the energy dependence?
A: Well, you got three points and two parameters for a straight line (in log-log
space), easy as hell right?  Report the dependence for each region too, may be
interesting spatial variation.  In principle can call Satoru if needed though
we want to work independently (note that he's in Japan, will be off phase ~12
hrs).

Next steps
----------
Profile fitting!  Rob around all week.
Spectra, when Brian can show how to merge the 750 ks of observations.
