Bugs observed at large B0, eta2; mu = 0
=======================================
Date: 2014 August 25

Whether this exact error occurs at other mu values is TBD...

(addendum, 2014 August 29 -- these errors are being identified by `models.py`,
mind you; the underlying code `FullEfflength_port.py` isn't throwing errors)

Initial observations
--------------------

Filament 1, finding upper bound on B0, using stderr causes massive overstep.
And, in general, the fitting (esp. at large B0) value just makes eta2 blow up
too fast.  Well, okay, no surprises.

But this resolution error I'm seeing below is really unexpected.

    Function call with B0 = 2203.167 muG; eta2 = 3187.787; mu = 0.000;
    rminarc = [ 60.  60.  60.]
		[ 1.90655254  1.61373412  1.20575015]
	Function call with B0 = 2203.167 muG; eta2 = 4220.534; mu = 0.000;
	rminarc = [ 60.  60.  60.]
		Resolution error! at [ 2.]
		[ 1.91005909  1.61616922  0.        ]

YES, this is reproducible and has some systematic behavior.  As I ramp up eta2,
more and more energy bands fall prey to the bug.  I guess that's a good sign,
that something is amiss and should be caught.

    >>> import numpy as np
    >>> import models_all as models
    >>> import models_all_exec_rewrite as mae
    >>> import snr_catalog as snrcat

    >>> f = mae.Fitter(snrcat.make_SN1006(), np.array([0.7, 1.0, 2.0]),
                       np.array([35.5, 31.94, 25.34]),
                       np.array([1.73, .97, 1.71]), 'tab')

    >>> models.full_width(f.snr, f.kevs, 0.0, 3187.787, 2203.167e-6)
    [ 1.9065526   1.61373416  1.20575018]
    >>> models.full_width(f.snr, f.kevs, 0.0, 3187.787, 2203.167e-6, rminarc=3)
    [ 1.98160087  1.65946649  1.17435883]

OOF, errors are order 3-7% with default SN1006 resolutions!
The FWHMs are just too narrow for rminarc = 60.

    >>> models.full_width(f.snr, f.kevs, 0.0, 4220.534, 2203.167e-6, rminarc=3)
    [ 1.98451644  1.66150468  0.        ]

Resolution error!  This is really unexpected.  We expect that with larger eta2,
the FWHMs will be larger.  But the whole thing just disappears ?!
(!!!) this will prevent fits from converging properly, and prevent root finding
from working correctly.  Damnation.
Example (Filament 1, upper bound on B0, mu=0, 3 energy bands):

    Function call with B0 = 862.015 muG; eta2 = 1544.835; mu = 0.000;
    rminarc = [ 60.  60.  60.]
		[ 8.03680808  6.76010475  4.77354712]
	Function call with B0 = 862.015 muG; eta2 = 2428.126; mu = 0.000;
	rminarc = [ 60.  60.  60.]
		[ 8.0839789   6.79035685  4.78658061]
	Function call with B0 = 862.015 muG; eta2 = 2428.127; mu = 0.000;
	rminarc = [ 60.  60.  60.]
		[ 8.08397893  6.79035687  4.78658062]
	Function call with B0 = 862.015 muG; eta2 = 3511.341; mu = 0.000;
	rminarc = [ 60.  60.  60.]
		Resolution error! at [ 2.]
		[ 8.10423357  6.80516313  0.        ]

Here the code tries to increase eta2 to get better FWHMs, closer to data.
But, runs into this strange/annoying error.
Therefore it thinks that B0 = 862.015 muG gives a TERRIBLE fit, when in fact
there may be an okay fit in parameter space, for even larger eta2.

ANOTHER very strange occurrence:
    Function call with B0 = 324.985 muG; eta2 = 3298.089; mu = 0.000;
    rminarc = [ 60.  60.  60.]
		[ 35.17726269  29.43171963  20.9445493 ]
    Function call with B0 = 324.985 muG; eta2 = 3301.938; mu = 0.000;
    rminarc = [ 60.  60.  60.]
            Box error! at [ 2.]
            [  35.17750267   29.43188739  899.66139881]
		
(?!) Python full model code didn't print an error as it normally would have done
And, usually the full width error would be 900. right on the dot.

This bug (these bugs?) drive the rootfinder nuts.  What the hell is wrong?!
Okay, this has to be fixed.  Sigh.



