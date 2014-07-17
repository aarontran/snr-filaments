Monday 2014 July 14
===================

Some material addressed at all-hands meeting (with Hiroya Yamaguchi, now, for
Nina's project), and a short afternoon meeting with Brian.

FWHM calculations
-----------------
* Present all the calculations, errors following Sean's paper for now.
  Mention issues in email.

Model code
----------
As you're figuring out the code -- write up that workflow thing, send to Rob
and Brian -- and we can send around to verify our understanding of the
math/code.

Misc. questions
---------------
* On computing `m_E` point-to-point, why use lower energies, instead of, say,
  using median energy of band (median energy of photon population -- of
  course that will be pretty close to lower energy, I guess, esp. in the
  farther spectrum)
  Brian: doesn't matter -- remember that the synchrotron spectrum falls off
  FAST, a median value will be like close to the lower energy anyways...
* grppha -- how do you select bin size (just eyeball, guesstimate?  specextract
  looks to be doing something more sophisticated, trying to equalize bin widths
  or something.  Seems like multiple criteria we could use)
  Brian: doesn't really matter.  Binning at 15 is fine.
* merge\_obs -- what energy to do exposure map at, for the bands?
  Brian: typically Brian just uses the middle energy.  Telescope response
  doesn't fall off so fast, so whatever.

