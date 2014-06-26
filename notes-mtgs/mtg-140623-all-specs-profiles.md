Monday 2014 June 23 (9:30-11a)
==============================

* Questions for Satoru (Brian will send today)
    - Errors from intensity (erg/cm^2/s)
    - FWHM calculation (where to measure FWHM? bottom of upstream side?)
    (what about downstream side?) (Rob commenting that one might be more
    conservative?)
    - Any reason that interpolating splines / similar weren't used?

* Try the Ressler model fit, see what happens

* Cut data from 1.7 - 2 keV, for bands -- guarantee nonthermal emission
* Brian will send new mosaics

* Quantify the effect of Si line emission on fits (do this today)
* (just ignore 1.7-2 keV, see how fits improve; or, add a Gaussian)
* (note to self: would be nice to get list of parameters...)
* XSPEC fitting (how to set limits on parameters)
    newpar 2 1.85 0.001 1.8 1.9 1.7 2.0
    # 1.85 is value
    # 0.001 is change in value
    # 1.8, 1.9 are soft limits
    # 1.7, 2.0 are hard limits

[addendum from short drop-ins with Rob/Brian on Monday/Tuesday:
XSPEC has equivalent width calculation functionality, just use that!)
