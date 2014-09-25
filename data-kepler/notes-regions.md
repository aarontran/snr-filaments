Kepler region selections
========================

Regions-2
---------
First, evaluate the quality of regions-1 and cull/split appropriately.

### Main ESE filament
* 22 -- spectrum mildly contaminated... but, few counts.
        structure is present but iffy, curiously shaped/curved
        Verdict: CUT
* 1 -- CUT (contamination)
* 2 -- extend northwards/southwards (get some region 1 cts), split in half
* 3 -- shrink north/south slightly, give cts to regions 2/4

* 4 -- extend north/south, split in half
* 5 -- move south, don't need to resize
* 6 -- shrink north slightly

  All clean.  Goal is to average ~200 cts at peak in 2-7 keV, for each region
  along this ESE filament (Kepler's ear, as Brian calls it)

* 7 -- contaminated, CUT

For reference:
> Region 02: length=38.32px, thickness=16.00px
>   Band 0.7-1kev: peak cts 146.0 in fit domain
>   Band 1-2kev: peak cts 664.0 in fit domain
>   Band 2-7kev: peak cts 279.0 in fit domain
> Region 03: length=39.17px, thickness=13.00px
>   Band 0.7-1kev: peak cts 160.0 in fit domain
>   Band 1-2kev: peak cts 750.0 in fit domain
>   Band 2-7kev: peak cts 407.0 in fit domain
> Region 04: length=37.55px, thickness=10.00px
>   Band 0.7-1kev: peak cts 125.0 in fit domain
>   Band 1-2kev: peak cts 623.0 in fit domain
>   Band 2-7kev: peak cts 266.0 in fit domain
> Region 05: length=35.63px, thickness=7.00px
>   Band 0.7-1kev: peak cts 164.0 in fit domain
>   Band 1-2kev: peak cts 397.0 in fit domain
>   Band 2-7kev: peak cts 194.0 in fit domain
> Region 06: length=29.88px, thickness=9.00px
>   Band 0.7-1kev: peak cts 119.0 in fit domain
>   Band 1-2kev: peak cts 526.0 in fit domain
>   Band 2-7kev: peak cts 271.0 in fit domain
> Region 07: length=30.18px, thickness=8.00px
>   Band 0.7-1kev: peak cts 95.0 in fit domain
>   Band 1-2kev: peak cts 375.0 in fit domain
>   Band 2-7kev: peak cts 156.0 in fit domain


### SE arc filaments
* 8 -- pretty clean.  0.7-1keV looks really bad.  spectrum chi2red ~ 1.25
       I shaved off about ~1 pixel thickness on north/south side,
       but stayed closer to south -- from looking at profiles, the 0.7-1keV
       band has more rim-shape to south, vs. north peaks a bit behind -- maybe
       due to thermal knot.
* 9 -- pretty clean, okay. chi2red ~ 1.32, a bit too much...
       Same deal. Trimmed 1/2 pixel thickness on north/south.

Region 08: length=29.29px, thickness=8.00px
  Band 0.7-1kev: peak cts 81.0 in fit domain
  Band 1-2kev: peak cts 389.0 in fit domain
  Band 2-7kev: peak cts 209.0 in fit domain
Region 09: length=29.63px, thickness=15.00px
  Band 0.7-1kev: peak cts 97.0 in fit domain
  Band 1-2kev: peak cts 337.0 in fit domain
  Band 2-7kev: peak cts 175.0 in fit domain

* 10 -- kinda clean, okay.  But, 1-2 keV profile is really meh, and we
        have fewer counts to work with.  Try narrowing/adjusting.
        let's cut this (see region 12)
* 11 -- contaminated, CUT. might work if we can get narrower FWHM...
* 12 -- same problems as 10.  NO 1-2 keV peak to speak of, but okay spectrum
        for consistency/simplicity, let's cut this
        someone else can work out the detailed kinematics/structure at kepler's
        shock, one day, when subarcsecond optics + new x-ray calorimeters get
        into the sky.

### Southern filaments, edges
* 13 -- clean.  But, profile has a double filament... ugh.  It messes up the
        fit.
* 14 -- clean, though not amazing counts.

Threw out 13, swung 14 around to get to thinner filament shape (away from
amorphous blob)

Region 13: length=35.17px, thickness=16.00px
  Band 0.7-1kev: peak cts 98.0 in fit domain
  Band 1-2kev: peak cts 338.0 in fit domain
  Band 2-7kev: peak cts 179.0 in fit domain
Region 14: length=32.54px, thickness=14.00px
  Band 0.7-1kev: peak cts 64.0 in fit domain
  Band 1-2kev: peak cts 270.0 in fit domain
  Band 2-7kev: peak cts 103.0 in fit domain

* 15 -- CUT (this inside arc is bright, but all contaminated.  too bad)
* 16 -- CUT
* 17 -- CUT

* 18 -- CUT
* 19 -- CUT

### Northern bits/pieces
* 20 -- eh. contaminated.  Try narrowing, avoid thermal stuff to right side.
        If that doesn't work, cut.  Might cut later for consistency, anyways.
        Good thing is that this one has more counts.
        Check the fits -- skimming Vink's paper on kepler kinematics, suggested
        northern regions might be dominated by thermal bremsstrahlung.

Region 20: length=33.29px, thickness=16.00px
  Band 0.7-1kev: peak cts 123.0 in fit domain
  Band 1-2kev: peak cts 319.0 in fit domain
  Band 2-7kev: peak cts 137.0 in fit domain

* 21 -- iffy.  faint structure, few counts. fair contamination, esp. below the
        1 keV mark.  I think I'll cut this.

### Results

11 regions selected.  7 slice up Kepler's ESE ear, 2 on another southern
filament, and 2 random ones at random spots


