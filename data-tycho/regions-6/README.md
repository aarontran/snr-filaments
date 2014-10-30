README
======

Tycho regions-6

Almost equivalent to Tycho regions-5, but with 20 instead of 22 regions
(21, 22 cut) and FWHMs manually culled/blacklisted where data are too poor

None of 0.7-1keV FWHMs for Regions 3, 5, 11, 19 are accepted if we impose a
data cap (i.e., FWHM must be measurable at 1/2 of data max).
So we are justified in throwing these guys out.
Furthermore, this step doesn't change our results substantially.
It is a bit of an eyeball thing, but it will do.  We see similar scatter in
many other bands, anyways.

Region peak information from `prep_profile_fit.py`:

    $ python ~/snr-research/code/profiles/prep_profile_fit.py \
        ../regions-6.physreg prfs prf-cts prf-proc \
        -l 0.7-1kev 1-1.7kev 2-3kev 3-4.5kev 4.5-7kev \
        -b 1 -w hanning -n 21 -v

    Region 01: length=51.09px, thickness=47.00px
      Band 0.7-1kev: peak cts 87.0 in fit domain
      Band 1-1.7kev: peak cts 597.0 in fit domain
      Band 2-3kev: peak cts 246.0 in fit domain
      Band 3-4.5kev: peak cts 152.0 in fit domain
      Band 4.5-7kev: peak cts 67.0 in fit domain
    Region 02: length=36.95px, thickness=17.00px
      Band 0.7-1kev: peak cts 79.0 in fit domain
      Band 1-1.7kev: peak cts 554.0 in fit domain
      Band 2-3kev: peak cts 179.0 in fit domain
      Band 3-4.5kev: peak cts 143.0 in fit domain
      Band 4.5-7kev: peak cts 48.0 in fit domain
    Region 03: length=37.24px, thickness=14.00px
      Band 0.7-1kev: peak cts 70.0 in fit domain
      Band 1-1.7kev: peak cts 462.0 in fit domain
      Band 2-3kev: peak cts 177.0 in fit domain
      Band 3-4.5kev: peak cts 113.0 in fit domain
      Band 4.5-7kev: peak cts 69.0 in fit domain
    Region 04: length=54.24px, thickness=24.00px
      Band 0.7-1kev: peak cts 132.0 in fit domain
      Band 1-1.7kev: peak cts 765.0 in fit domain
      Band 2-3kev: peak cts 307.0 in fit domain
      Band 3-4.5kev: peak cts 199.0 in fit domain
      Band 4.5-7kev: peak cts 89.0 in fit domain
    Region 05: length=48.90px, thickness=24.00px
      Band 0.7-1kev: peak cts 115.0 in fit domain
      Band 1-1.7kev: peak cts 872.0 in fit domain
      Band 2-3kev: peak cts 331.0 in fit domain
      Band 3-4.5kev: peak cts 211.0 in fit domain
      Band 4.5-7kev: peak cts 93.0 in fit domain
    Region 06: length=33.58px, thickness=26.00px
      Band 0.7-1kev: peak cts 140.0 in fit domain
      Band 1-1.7kev: peak cts 919.0 in fit domain
      Band 2-3kev: peak cts 319.0 in fit domain
      Band 3-4.5kev: peak cts 228.0 in fit domain
      Band 4.5-7kev: peak cts 86.0 in fit domain
    Region 07: length=36.30px, thickness=22.00px
      Band 0.7-1kev: peak cts 107.0 in fit domain
      Band 1-1.7kev: peak cts 811.0 in fit domain
      Band 2-3kev: peak cts 322.0 in fit domain
      Band 3-4.5kev: peak cts 207.0 in fit domain
      Band 4.5-7kev: peak cts 76.0 in fit domain
    Region 08: length=44.67px, thickness=13.00px
      Band 0.7-1kev: peak cts 80.0 in fit domain
      Band 1-1.7kev: peak cts 580.0 in fit domain
      Band 2-3kev: peak cts 231.0 in fit domain
      Band 3-4.5kev: peak cts 189.0 in fit domain
      Band 4.5-7kev: peak cts 65.0 in fit domain
    Region 09: length=58.20px, thickness=13.00px
      Band 0.7-1kev: peak cts 76.0 in fit domain
      Band 1-1.7kev: peak cts 473.0 in fit domain
      Band 2-3kev: peak cts 219.0 in fit domain
      Band 3-4.5kev: peak cts 135.0 in fit domain
      Band 4.5-7kev: peak cts 69.0 in fit domain
    Region 10: length=58.48px, thickness=18.00px
      Band 0.7-1kev: peak cts 84.0 in fit domain
      Band 1-1.7kev: peak cts 482.0 in fit domain
      Band 2-3kev: peak cts 194.0 in fit domain
      Band 3-4.5kev: peak cts 161.0 in fit domain
      Band 4.5-7kev: peak cts 59.0 in fit domain
    Region 11: length=37.63px, thickness=12.00px
      Band 0.7-1kev: peak cts 48.0 in fit domain
      Band 1-1.7kev: peak cts 244.0 in fit domain
      Band 2-3kev: peak cts 107.0 in fit domain
      Band 3-4.5kev: peak cts 80.0 in fit domain
      Band 4.5-7kev: peak cts 32.0 in fit domain
    Region 12: length=37.50px, thickness=11.00px
      Band 0.7-1kev: peak cts 47.0 in fit domain
      Band 1-1.7kev: peak cts 290.0 in fit domain
      Band 2-3kev: peak cts 124.0 in fit domain
      Band 3-4.5kev: peak cts 110.0 in fit domain
      Band 4.5-7kev: peak cts 49.0 in fit domain
    Region 13: length=27.60px, thickness=27.00px
      Band 0.7-1kev: peak cts 131.0 in fit domain
      Band 1-1.7kev: peak cts 833.0 in fit domain
      Band 2-3kev: peak cts 318.0 in fit domain
      Band 3-4.5kev: peak cts 242.0 in fit domain
      Band 4.5-7kev: peak cts 104.0 in fit domain
    Region 14: length=51.36px, thickness=11.00px
      Band 0.7-1kev: peak cts 88.0 in fit domain
      Band 1-1.7kev: peak cts 483.0 in fit domain
      Band 2-3kev: peak cts 164.0 in fit domain
      Band 3-4.5kev: peak cts 107.0 in fit domain
      Band 4.5-7kev: peak cts 41.0 in fit domain
    Region 15: length=51.24px, thickness=11.00px
      Band 0.7-1kev: peak cts 130.0 in fit domain
      Band 1-1.7kev: peak cts 715.0 in fit domain
      Band 2-3kev: peak cts 285.0 in fit domain
      Band 3-4.5kev: peak cts 136.0 in fit domain
      Band 4.5-7kev: peak cts 63.0 in fit domain
    Region 16: length=44.58px, thickness=12.00px
      Band 0.7-1kev: peak cts 125.0 in fit domain
      Band 1-1.7kev: peak cts 800.0 in fit domain
      Band 2-3kev: peak cts 293.0 in fit domain
      Band 3-4.5kev: peak cts 189.0 in fit domain
      Band 4.5-7kev: peak cts 78.0 in fit domain
    Region 17: length=40.85px, thickness=12.00px
      Band 0.7-1kev: peak cts 139.0 in fit domain
      Band 1-1.7kev: peak cts 777.0 in fit domain
      Band 2-3kev: peak cts 291.0 in fit domain
      Band 3-4.5kev: peak cts 188.0 in fit domain
      Band 4.5-7kev: peak cts 62.0 in fit domain
    Region 18: length=48.71px, thickness=47.00px
      Band 0.7-1kev: peak cts 129.0 in fit domain
      Band 1-1.7kev: peak cts 691.0 in fit domain
      Band 2-3kev: peak cts 238.0 in fit domain
      Band 3-4.5kev: peak cts 124.0 in fit domain
      Band 4.5-7kev: peak cts 57.0 in fit domain
    Region 19: length=48.26px, thickness=19.00px
      Band 0.7-1kev: peak cts 68.0 in fit domain
      Band 1-1.7kev: peak cts 371.0 in fit domain
      Band 2-3kev: peak cts 134.0 in fit domain
      Band 3-4.5kev: peak cts 104.0 in fit domain
      Band 4.5-7kev: peak cts 39.0 in fit domain
    Region 20: length=44.03px, thickness=14.00px
      Band 0.7-1kev: peak cts 80.0 in fit domain
      Band 1-1.7kev: peak cts 454.0 in fit domain
      Band 2-3kev: peak cts 195.0 in fit domain
      Band 3-4.5kev: peak cts 132.0 in fit domain
      Band 4.5-7kev: peak cts 44.0 in fit domain


