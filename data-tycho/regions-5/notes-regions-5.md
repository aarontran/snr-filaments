Notes on picking regions 5
==========================
2014 September 22

Goal: split up the regions and get more FWHMs in each filament
Look at the profiles (inspect images), look at the max number of counts for
each region, and determine which regions I can slice up.
All of these regions/filaments have already had their spectra well-vetted.

Now is the time to look back at `regions-all.reg`... will do so more thoroughly
on next round.  Right now, goal is to try getting more profiles and more fit
numbers.

Numbers printed below are for regions-4 (old numbering system).  I changed up
the numbering system with regions-5

Region 01: length=61.46px, thickness=29.00px
  Band 0.7-1kev: peak cts 132.0 in fit domain
  Band 1-1.7kev: peak cts 849.0 in fit domain
  Band 2-3kev: peak cts 385.0 in fit domain
  Band 3-4.5kev: peak cts 270.0 in fit domain
  Band 4.5-7kev: peak cts 111.0 in fit domain

  Excepting 0.7-1 keV band (which is still pretty clean),
  go ahead and try slicing this up.
  Shrink Region 13 just slightly, to make more space.

Region 02: length=31.25px, thickness=22.00px
  Band 0.7-1kev: peak cts 92.0 in fit domain    (blacklisted)
  Band 1-1.7kev: peak cts 507.0 in fit domain
  Band 2-3kev: peak cts 223.0 in fit domain
  Band 3-4.5kev: peak cts 172.0 in fit domain
  Band 4.5-7kev: peak cts 76.0 in fit domain

  We are running low on counts in 4.5-7 keV.
  Extend a little northwards, make it a bit longer.
  THEN, try to slice it into two -- aim for proportional numbers of cts.

Region 03: length=27.60px, thickness=27.00px
  Band 0.7-1kev: peak cts 131.0 in fit domain
  Band 1-1.7kev: peak cts 833.0 in fit domain
  Band 2-3kev: peak cts 318.0 in fit domain
  Band 3-4.5kev: peak cts 242.0 in fit domain
  Band 4.5-7kev: peak cts 104.0 in fit domain

  0.7-1 keV is a little iffy, as is 4.5-7 keV.
  Not much room to extend, as w/ Region 2

Region 04: length=51.60px, thickness=20.00px
  Band 0.7-1kev: peak cts 180.0 in fit domain
  Band 1-1.7kev: peak cts 1053.0 in fit domain
  Band 2-3kev: peak cts 379.0 in fit domain
  Band 3-4.5kev: peak cts 231.0 in fit domain
  Band 4.5-7kev: peak cts 93.0 in fit domain

  Split, but extend south slightly -- make sure we have enough counts at 4.5-7
  keV

Region 05: length=43.89px, thickness=25.00px
  Band 0.7-1kev: peak cts 271.0 in fit domain
  Band 1-1.7kev: peak cts 1609.0 in fit domain
  Band 2-3kev: peak cts 588.0 in fit domain
  Band 3-4.5kev: peak cts 413.0 in fit domain
  Band 4.5-7kev: peak cts 152.0 in fit domain

  YES, split this in half!  Don't extend.  Plenty of counts.

Region 06: length=40.71px, thickness=47.00px
  Band 0.7-1kev: peak cts 132.0 in fit domain   (blacklisted)
  Band 1-1.7kev: peak cts 691.0 in fit domain
  Band 2-3kev: peak cts 238.0 in fit domain
  Band 3-4.5kev: peak cts 121.0 in fit domain
  Band 4.5-7kev: peak cts 57.0 in fit domain

  4.5-7 keV is at limit of rattiness.  This one is already VERY faint, and
  spectrum fit is not great.

Region 07: length=40.31px, thickness=26.00px
  Band 0.7-1kev: peak cts 123.0 in fit domain
  Band 1-1.7kev: peak cts 709.0 in fit domain
  Band 2-3kev: peak cts 269.0 in fit domain
  Band 3-4.5kev: peak cts 203.0 in fit domain
  Band 4.5-7kev: peak cts 74.0 in fit domain

  This one is hard -- sitting right on a lump of filaments!
  But, try extending left (touch Region 6) and splitting.

  ALSO -- in accordance w/ adding regions where 0.7-1 keV is shot,
  add at least one more region farther southeast -- at 3 filament splay

Region 08: length=51.09px, thickness=47.00px
  Band 0.7-1kev: peak cts 87.0 in fit domain
  Band 1-1.7kev: peak cts 597.0 in fit domain
  Band 2-3kev: peak cts 246.0 in fit domain
  Band 3-4.5kev: peak cts 152.0 in fit domain
  Band 4.5-7kev: peak cts 67.0 in fit domain

  At rattiness limit, and FWHMs very wide.
  Spectral fit is not great (redchisqr ~ 1.2).

  ALSO -- try adding back the region in between, and possibly look at wispier
  filaments to south (bad in 0.7-1 keV for sure).
  (NOPE -- the region in between ain't gonna work.  Too smeared out in 1-2 kev)

Region 09: length=37.29px, thickness=32.00px
  Band 0.7-1kev: peak cts 153.0 in fit domain
  Band 1-1.7kev: peak cts 991.0 in fit domain
  Band 2-3kev: peak cts 374.0 in fit domain
  Band 3-4.5kev: peak cts 279.0 in fit domain
  Band 4.5-7kev: peak cts 115.0 in fit domain

  Extend slightly more south and split.  Enough counts for sure.
  Correction -- nope.  This runs into the problem of smearing out 1-2 kev.
  Well, try splitting anyways and see what happens.

Region 10: length=54.77px, thickness=31.00px
  Band 0.7-1kev: peak cts 164.0 in fit domain
  Band 1-1.7kev: peak cts 1033.0 in fit domain
  Band 2-3kev: peak cts 415.0 in fit domain
  Band 3-4.5kev: peak cts 246.0 in fit domain
  Band 4.5-7kev: peak cts 104.0 in fit domain
Region 11: length=40.02px, thickness=30.00px
  Band 0.7-1kev: peak cts 158.0 in fit domain   (blacklisted)
  Band 1-1.7kev: peak cts 994.0 in fit domain
  Band 2-3kev: peak cts 391.0 in fit domain
  Band 3-4.5kev: peak cts 270.0 in fit domain
  Band 4.5-7kev: peak cts 105.0 in fit domain
Region 12: length=37.72px, thickness=33.00px
  Band 0.7-1kev: peak cts 169.0 in fit domain
  Band 1-1.7kev: peak cts 1225.0 in fit domain
  Band 2-3kev: peak cts 453.0 in fit domain
  Band 3-4.5kev: peak cts 301.0 in fit domain
  Band 4.5-7kev: peak cts 131.0 in fit domain

  Taking 10,11,12 together -- carve 4 regions out of 3, and shrink Region 12
  since it has a few more counts than its neighbors.
  Aiming for something like 80 cts/region in 4.5-7 keV, and 120 cts/region in
  0.7-1 keV (at peak).

Region 13: length=44.67px, thickness=16.00px
  Band 0.7-1kev: peak cts 104.0 in fit domain
  Band 1-1.7kev: peak cts 694.0 in fit domain
  Band 2-3kev: peak cts 279.0 in fit domain
  Band 3-4.5kev: peak cts 237.0 in fit domain
  Band 4.5-7kev: peak cts 80.0 in fit domain

  Kind of pushing the rattiness limit.
  But, can afford to shrink just a tiny bit.
  Don't split.


