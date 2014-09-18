Bugs observed at large B0, eta2
=======================================
Date: 2014 August 25


These errors are being identified by `models.py`, mind you;
the underlying code `FullEfflength_port.py` isn't throwing errors)

I've observed this at mu = 0.0, 0.5, so it's not limited to mu = 0.
I've no idea if it appears at mu=1 or mu>1, but my guess is that the fits
simply don't explore that area of parameter space at large mu.

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


Debugging log
-------------
2014 August 30

I see two different types of errors popping up. From an attempt at error before,
it seems like there may be a transition between error types.
From `20140829_error_rootfinding.md`:

* B0 = 795.590 muG; eta2 = 146434.006; mu = 0.500 -- max eta2 that computes okay
* B0 = 795.590 muG; eta2 = 146434.216; mu = 0.500 -- box error at [ 2.]
* B0 = 795.590 muG; eta2 = 146893.065; mu = 0.500 -- box error at [ 2.]
* B0 = 795.590 muG; eta2 = 163307.996; mu = 0.500 -- resolution error at [ 2.]

So, I'll start working from this.  Disable adaptive rminarc, for now.
(it only makes it more complex)

### Boilerplate code

`
import numpy as np
import models
import models_exec as mex
import snr_catalog as snrcat
f = mex.Fitter(snrcat.make_SN1006(), np.array([0.7,1.,2.]), 1,1,1)
def fd(mu=0.5, eta2=146434, B0=795.590e-6, **kwargs):
    return f.width_full(mu, eta2, B0, irad_adapt=False, **kwargs)
%load_ext autoreload
%autoreload 2
`

### Useful inserted debug code

In case I need to debug FWHM calculation further in the future...

    rmin = spopt.bisect(f_thrsh, rmin_a, rmin_b)
    #rmin, rmin_res = spopt.brentq(f_thrsh, rmin_a, rmin_b, full_output=True)
    rmax = spopt.bisect(f_thrsh, rmax_a, rmax_b)
    #rmax, rmax_res = spopt.brentq(f_thrsh, rmax_a, rmax_b, full_output=True)

    #print 'DEBUG: finding half-max crossings'
    #print 'crossing indices: {}, {}'.format(inds_rmin, inds_rmax)
    #print 'min: search between {}, {}'.format(rmin_a, rmin_b)
    #print 'max: search between {}, {}'.format(rmax_a, rmax_b)
    #print 'min: intensity bracket {}, {}'.format(f_int(rmin_a), f_int(rmin_b))
    #print 'max: intensity bracket {}, {}'.format(f_int(rmax_a), f_int(rmax_b))
    #print 'halfpk intensity {}'.format(halfpk)
    #print ' min intensity gap {}, slope {}'.format(f_int(rmin_b)-f_int(rmin_a),
    #            (f_int(rmin_b)-f_int(rmin_a)) / (rmin_b - rmin_a))
    #print ' max intensity gap {}, slope {}'.format(f_int(rmax_b)-f_int(rmax_a),
    #            (f_int(rmax_b)-f_int(rmax_a)) / (rmax_b - rmax_a))
    #print ' brentq min debug: conv={}, fcalls={}, iter={}'.format(rmin_res.converged,
    #        rmin_res.function_calls, rmin_res.iterations)
    #print ' brentq max debug: conv={}, fcalls={}, iter={}'.format(rmax_res.converged,
    #        rmax_res.function_calls, rmax_res.iterations)
    #print 'rmin, rmax = {}, {}'.format(rmin, rmax)
    #print 'crossing intensities: {}, {}'.format(f_int(rmin), f_int(rmax))
    #print 'pk intensity {}'.format(pk)



### Debugging

In the case of box error, where we expect FWHM ~ 30 arcsec,
we find the crossings in `irmax` grid correctly.
And, we give brentq reasonable bracketing values.
But, rmin comes out as 0.0, rather inexplicably!

It seems to break between eta2 = 146434.2819*978*, 146434.281*979*

* eta2 = 146434.2819978
    DEBUG: finding half-max crossings
    crossing indices: [242], [397]
    min: search between 0.973666666667, 0.973833333333
    min: intensity bracket 1.59593490018e-159, 1.60281266752e-159
    halfpk intensity 1.59631509252e-159
     min intensity gap 6.87776733967e-162, slope 4.1266604038e-158
     brentq min debug: conv=True, fcalls=30, iter=29
    rmin, rmax = 0.973833333332, 0.999567926696
    crossing intensities: 1.60281266747e-159, 1.5459516023e-159

* eta2 = 146434.2819979
    DEBUG: finding half-max crossings
    crossing indices: [242], [397]
    min: search between 0.973666666667, 0.973833333333
    min: intensity bracket 1.59593489994e-159, 1.60281266728e-159
    halfpk intensity 1.59631509227e-159
     min intensity gap 6.87776733862e-162, slope 4.12666040317e-158
     brentq min debug: conv=True, fcalls=2, iter=1
    rmin, rmax = 0.0, 0.999567926696
    crossing intensities: 3.12711776556e-160, 1.54595160206e-159

Wow, I managed to break scipy's brentq.
The solution appears to be as simple as, use `scipy.optimize.bisect`.
Bisection is more robust, and for us doesn't drop performance because
our precomputed grid has already narrowed the search space greatly.

A little fiddling suggests that it may be the intensity issue.
The other test cases below fail at ~1e-200, 1e-170, 1e-160.
If I artificially scale up intensity numbers by 1e100, the errors disappear.
HMMMMMMM.

### Some test cases to check (copy-pasted from above)

* `B0 = 324.985 muG; eta2 = 3301.938; mu = 0.000; rminarc = [ 60.  60.  60.]`
    Box error! at [ 2.]
    [  35.17750267   29.43188739  899.66139881]

* `B0 = 862.015 muG; eta2 = 3511.341; mu = 0.000; rminarc = [ 60.  60.  60.]`
	Resolution error! at [ 2.]
	[ 8.10423357  6.80516313  0.        ]

* `B0 = 2203.167 muG; eta2 = 4220.534; mu = 0.000; rminarc = [ 60.  60.  60.]`
	Resolution error! at [ 2.]
	[ 1.91005909  1.61616922  0.        ]


Solution and explanation
------------------------

Replace brentq with bisect in my FWHM code.
Additionally, normalize gridded intensity and intensity function in my FWHM
computation.  Any change in results due to this should only be due to
1. differences in brentq/bisect tolerance, or
2. floating point error (from intensity rescaling?!)


Conclusions, future work
------------------------

Today (this week, rather), `scipy.optimize.brentq` failed us.
We put `brentq` to the test and, facing witheringly minuscule input,
`brentq` died in the line of duty.

I don't fully understand where or why it fails.  The most possibly relevant
numbers are:
* intensity ~ 1.60e-159
* intensity gap (between bracketing pts) ~ 6.88e-162
* slope ~ 4.127e-158

The other cases all have smaller values of intensity/gap/slope.  Perhaps the
relative change (change in intensity / intensity) matters.


Truth be told I'm not inclined to dig farther now and understand why.
To explore this further: roll back to this commit, remove my scalings and use
brentq again, and use these test cases.


Addendum: floating point limits
-------------------------------
Ramping up eta2 to 5e5 or so, the full model code is balking at me.
(eta2 = 4.3e5, B0 = 795.590e-6, mu = 0.5)

Intensity values are starting to hit the smallest representable numbers in
Python, around 1e-300 or so.  The code starts behaving weirdly, and eventually
you just get nans/infs/zeros out or whatever...


This might be solved by, e.g., removing some of the normalization in the
intensity calculation (or throwing on some arbitrary prefactor, somewhere).
But... these cases are so extreme, and may not even be physically meaningful,
sensible, or relevant. So I won't touch it.

But it did reveal that, there's no need to (un)normalize intensity in FWHM
calculations, it just makes the situation worse.

See also, in re smallest float: http://stackoverflow.com/a/6163157

