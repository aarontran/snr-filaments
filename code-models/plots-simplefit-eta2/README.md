README for SN1006 plots
=======================
(generated 2014 July 24)

Using SN 1006 data.

Procedure: grid over eta2 and mu.
At each point, use the simple model to compute best fit value of B0
(this traces out a line in eta2-B0 space, capturing best FWHMs)

Then, plot chi-squared as a function of eta2.  This tells us where the best fit
pair (eta2, B0) is for a given value of mu.

Do this for all mu values.  This shows us whether a particular mu value is
favored, and potentially why (but the physical interpretation is often unclear,
except in cases where diffusion is unimportant (eta2 small, mu no longer
relevant)).



Purpose of this is to get a handle on parameter behavior, and how model gives
reasonable FWHM values.

At fixed eta2 and mu, fitting for B0 is reasonably robust.

But, it is more difficult to identify the best fit eta2.

And, the best fit mu is probably not even significant / useful with so few data
points.



This was generated with the 2014 July 24 version of `models_all.py` (SN1006
values, widths + errors taken from original `Widthfun.py` from SEan).

Code is being modified to focus on tabulating values for multiple B0, and not
generate fits.  But to reproduce this is not so hard.
