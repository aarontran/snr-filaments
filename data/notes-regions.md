Notes on handpicked regions for Tycho SNR profiles
==================================================
(last modified: 2014 Jun 26)

This is generally an append-only log, to record what was done for region
selection (what regions were added/discarded and on what grounds).

It's messy but should give a general idea.  In general, work with
regions-all, then cut that down to a given set of regions (version 1,2,3,4,
whatever).

Procedure, notes, output
========================
Start by issuing shell command:

    ds9 -rgb \
        -red 0.7-1kev_mosaic.fits \
            -scale limits 4e-8 2.1e-6 \
        -green 1-2kev_mosaic.fits \
            -scale limits 4e-8 2.42e-6 \
        -blue 2-7kev_mosaic.fits \
            -scale limits 3e-8 1.5e-6 \
        -rgb lock scale yes \
        -histequ

Rough procedure for eyeballing regions:
* Load using the command above.
* Currently using scale=hist and looking at 2-7 keV profiles
* Draw from inside out, so plots increase in radial distance
* Eyeball, e.g. by lining up edge of box with rim (R. Petre)
* Make note of bad/contaminated regions too, maybe some way to remove
  thermal contamination?
* Save regularly!

In practice, I find myself relying mostly on maximizing peak height.  Easiest
to do by eye.  Also trying to minimize contamination, try to minimize the
trough behind the peak (between peak and thermal emission tufts)
I have not been thinking about projection when choosing filaments

When thinking about region sizes/counts:
* SN 1006 has diameter 0.5 deg., vs. ~0.14 deg. for SN 1572
* But count rates are much higher for SN 1572

Output: ds9 .reg file (can also save to CIAO formats or other things)
Saved ds9 region files are plaintext and contain width information
Direction of "thickness" is determined by ordering of two stored points,
so region's shape and orientation are fully specified.

Region numbering and comments
=============================

Starting from the westernmost (lowest RA, rightmost in ds9) rim of the shell,   
I move counterclockwise about the shell, hand-selecting regions and making
notes.  Regions are numbered in order.

[good/blah/bad!], n) gives an idea of (1) region quality, (2) number of
filaments.

The regions are color-coded by quality to assist in follow-up analysis.

*Northwest corner*

Start RA,dec = (00:24:40.6, +64:07:52).

01. (good, 1) Region chosen to fit between two little blue knots
02. (bad!, 2) Fit between two little knots.  Looks like two faint filaments.

03. (good, 1) Pretty clear, some contamination
04. (blah, 1) Strong signal, but moderate contamination in back
05. (good, 1) Pretty clear, some contamination

Regions 3-5 are the same filament, but 4 captures middle contamination.
Originally a single region, which was too smeared out.  Three split regions are
all better than original region.

06. (blah, 1) Looks like some other filaments sneak in, towards the north
              But I stay away from those, and try to minimize contamination
              (again, avoid a small knot of emission towards the north)
07. (bad!, 1) Thin, weak filament, bad contamination.  Stronger in middle,
              but more contamination in middle too.
08. (blah, 2) NW filaments separated into two peaks
09. (good, 2) Bright, interacting NW filaments, far away from thermal emission.
              Treat as one strong filament? (read abt modeling process)
              The region has some mush to the south, as dim filament emerges
              but, it's not so strong. I tried to focus on the bright filament
10. (bad!, 2-3) Interacting filaments, above bright NW.
                North, single filament, but contamination is strong
                South, less contamination, but two(three?) filaments

*Northeast corner*

11. (good, 1-2) NE rim filaments, looks like a merger of two filaments on
                southern side (counterclockwise side).  Bright edge.
                Weak signal, but VERY low contamination.
12. (bad!, 1-2) South of nice NE rim, where the filaments begin splaying apart
                A very narrow selection, where the filaments split but are
                still very close together.  Not so good.
13. (bad!, 1-2) Further along NE rim, where one filament has split out far
                enough so we can see two distinct peaks (+ thermal at edge)
                Also very narrow selection.
14. (blah, 1-2) Further along NE rim, Two overlapping filaments.
                Previous split filament has merged with thermal radiation.
                Clean, but unfortunately kind of broad due to two filaments.
15. (good, 2) Clear-ish front filament, with thermal emission in back +
              a splayed off filament.  Tried to select region where filaments
              were reasonably distinct.
16. (bad!, 1-2) Front filament is just being hit with contamination.
                Back (2nd filament) is overrun with thermal stuff, but is still
                visible!  Maybe useful for tests to remove contamination.

In northeast corner (between regions 16, 17), there's a strong
point-like source in 2-7kev band (strongest src in histogram, decently strong
in other bands too).  Avoid this area.

17. (bad!, 2-?) Multiple filaments!  Speculating, perhaps multiple interactions
                with gas cloud?
                I think I can see up to 4 filaments in the profile.
                At the back is a tuft of thermal emission.
18. (blah, 1?) Looks like one nice filament, but contaminated.
               Tuft of green (thermal?) stuff will contaminate local background
               estimation.  So I extended region out a bit farther.
               If we use this, address local background on this one specially.
19. (good, 1-2) Strong signal!
                I have extended it enough to pick up a 2nd filament.
                Although I could exclude it -- it's so strong, that it becomes
                clearer just by adding its many counts, even as we pick up
                additional stuff on the tail.
                Best viewed with sqrt/asinh scaling, mush on histequ
20. (blah, 1-2) Strong, but the two filaments are too close together
                So they mush up, and there's stuff on the tail...

*Southeast corner*

21. (good, 1?) Clear filament, with slight contamination.
               Not especially strong (esp. compared to 19/20).
               Shock does smear out toward to south, so I make the region a
               little narrower.
22. (bad!, 1) Heavily contaminated filament, and very weak to begin with.
23. (bad!, 1) Same filament as in 22, but fainter and with less contamination.
              Arguably could be two filaments... but really can't say.
              These two (22, 23) are the very faint SE filament, which is best
              viewed in sqrt scaling (hard to see in asinh or histequ).
24. (blah, 2) Two relatively faint, but not too badly contaminated filaments.
25. (blah, 2) Two overlapping filaments, I tried to get just their merger.
              Okay but the peak is a bit broader with two.
              Some thermal emission contaminates, but not too bad.
26. (blah, 2) Two overlapping filaments, a bit stronger, but also more thermal
              contamination present.  Not great.

*Southwest corner*

27. (bad!, 1) Very faint, smeared out shock.  Doesn't appear too contaminated,
              but I cannot tell for sure.  Visible front peak in Green/Blue,
              less visible or not present in red; red thermal emission fluffs
              up in the back.
28. (blah, 2) Same problems as region 27, but stronger signal and peaks appear
              more separated.  Narrowed integration region to get stronger
              signal.
29. (bad!, 1-2) Decent front peak, but then followed by mush... does not appear
                too contaminated though, just smearing in the back from a
                possible second filament.
30. (bad!, 1?) Contaminated in the back, kind of weak signal.  But, definitely
               a filament.
31. (bad!, 1) Large blob in filament, with thermal emission not far behind.
              Pretty contaminated.  Peak is fairly wide, not sure why.
              Peak width could be meaningful?
              Or (e.g. Chevalier 1977, Figs 1, 2), it could be a blob of ejecta
              linked to the Rayleigh-Taylor instability shape right behind it?
32. (good, 1) Long SW filament, fairly uncontaminated.  Sharp rise in the back
              from thermal stuff, but we can ignore it.  The key features of
              the filament are unaffected.
33. (blah, 1) A bit more contamination.  The peak is a little broader,
              the back of the peak more difficult to make out.
34. (good, 1) Slight contamination (similar to 33), but peak is much stronger
              thus we are saved.  Very sharp peak.

Regions 32-34 are all the same filament, but broken into parts.
Motivated by (1) shape of thermal emission behind shock, (2) varying strength
of filament emission.

35. (blah, 1) Filament blob -- some kind of (1) spreading, or
              (2) ejecta interacting with shock wave.
              Contamination pretty close by.
              Getting this region for completeness, esp. if spreading is
              physically significant.
              It looks as if the shock wave changes direction here, or as if
              there are two filaments interacting here.

36. (good, 1) Strong peak, low contamination
37. (good, 1) Strong peak, low contamination but worse than in 34
              (contamination is closer to filament)
38. (good, 1) Strong peak, very low contamination!
39. (blah, 1) Good peak, but moderate contamination (close to filament)
40. (good, 1) Good peak, lower contamination (not as good as region 1)
  
Regions 36-40 and 1 are all part of the same, very bright filament.

*Some notes on gauging signal strength, contamination*

Generally, bad!/blah regions are contaminated in the back and/or have multiple,
overlapping filaments.  However, they all have visible front peaks.
In some places I comment about the peak's strength, which can sometimes rise
out strongly from the background contamination.

Qualitatively, very strongest signals are around 7e-7 count average (peak).
Troughs around 1e-7 to 5e-8 are good, depending on peak strength
Background outside SNR is usually ~1e-8
Example of a decent faint one (22), peak is ~8e-8, trough 4e-8.

Most regions are fainter/weaker in red, green (mainly in red), but most appear
to be usable still.


Region classification
=====================

After picking out my band of 40, I need to pare them down...

*First pass, using binned data*

Good single filaments: 1, 3, 5, 9, 11, 15?, 19, 21, 32,33,34, 36,37,38, 40

I'm now going through the good filaments only, to verify that they are indeed
good.

Accept: 1, 9, 11, 21, 32, 34, 36, 37, 38, 40
3, 4, 5 are uncertain -- some contamination, but visible peaks
20, 25 are strong enough that they could be workable...

Reject:
* 3, 4, 5 (contamination?)
* 15 (multiple filaments)
* 19, 20 (multiple filaments, contamination)

(1,36,37,38,40; 9; 11; 21; 32,34 are distinct... so only 5 distinct filaments)

I generate a new region file with the good ones only, for testing...

*Second pass, using unbinned data*

Now I've regenerated the profile plots, I go through one more time to ensure I
have good regions. I look at 1) unbinned profiles, 2) DS9 image, 3) fiddle with
projections in DS9 to double check what I'm seeing.

Idea: plots are good for cutting out if they don't look very fittable.
But, thermal contamination needs spectra to assess.

Rough criteria -- single peak, no obvious contamination in RGB image, red peak
should be clear down to ~half max.
By "good", I mean good for profile fitting extraction.  Some of the bad ones,
may be interesting to generate spectra for anyways (mainly southern rim).

GOOD: 1, 3, 5, 9, 11, 36-38, 40 (unequivocally good)
????: 32-33, 34, 35 (likely good, but peaks a bit fainter)
????: 12, 18, 19, 21 (iffy, contaminated)
BAD!: 2, 4, 6, 7, 8, 10, 13-17, 20, 22-24, 25-26, 27-31, 39 (don't use)

From spectra for some good regions:
1, 11 are slightly contaminated, fiddle those
9 is exceptionally clean (clear from image)
32 looks pretty clean, but noisy
34 is slightly contaminated
36-38, 40 are less contaminated than (1, 11, 34), but kind of noisy

*Qualitative summary*
Bit of WNW band (regions 3-5) is okay, though maybe contaminated
Bright NW band (region 9) is very good
Weak edge of NE band is good (region 11), but, flays apart too much farther
down (12-17)
ENE, ESE rims (regions 18-21) are strong but contaminated.  I'll try a few.
Wispy SE rims (regions 22-24) are hopeless
SSE rims (regions 25-26) are hopeless for profile fitting, very smeared in red
SSW rims (regions 27-31) have same problem, smushed out in red.
W filaments (regions 32-40, 1) are excellent, though some areas of smushing
and knots

*Improvements*

* Twiddle with regions 3-5, 11-12 ish.
* Region 18 -- make narrower, smaller, try again.
* Region 19/20 COULD be workable, must be right between the knots, where there
* are two filaments.  try again.
* Region 21, try again but without the stuff in the back (cutback and adjust)
* Region 32-33, try again, get a wider area, avoid thermals.
* Region 34, also twiddle.

* Regions 13-17 not useful, but try to get spectra (is it all synchrotron?)
* Regions 25-26 not useful, but try to get spectra -- strictly avoid thermals
* Same for 27-30 area
(could assess projection effects...)

Generate this in `profiles_all_2.reg` or something similarly named

Region selection, round 2
=========================

Start working on `profiles_all_2.reg`.
1. Throw away bad spectra -- we can return to `profiles_all.reg`.
2. Only the good/useful ones for now.  Southern rim, revisit later.
3. Apply cutbacks on all regions, extend backwards soon.
4. Look at RED band for profiles, most discriminatory

The single region 11, I have now divided up between 11 and 12.

Region in between 11-12, looks very bad in RED band -- split into two
filaments, lots of mush.

32,33 are just MUSH in red.  I let them overlap a bit, to test them out.

(edit [jun 26] -- these are derived from / modified from `regions-all`)

Region evaluation (of round 2 picks)
------------------------------------

Yes, the numbering is horrible.
For regions where I have not selected spectra, refer to profile plots for
`regions-all.reg`, which will illustrate rim mushiness or splaying.

Warning -- reduced-chi-squared is actually not that helpful.
Best way to tell spectrum quality is by looking for the silicon line.
Remember Brian's suggestions: nH 0.4 to 0.8, index 2 to 3, chi2 < ~1.5.

New numbering, old numbering, profile/RGB image quality, spectrum quality

1. (1) [GOOD]       (nH = 0.71, index = 2.86; chi2 = 1.05)  (small line)
2. (3) [okay?]      (nH = 0.67, index = 2.57; chi2 = 1.08)  (clean)
3. (5) [okay?]      (nH = 0.58, index = 2.69; chi2 = 1.08)  (clean)
4. (9) [GOOD]       (nH = 0.56, index = 2.80; chi2 = 1.14)  (AMAZING)
5. (11) [okay?]     (nH = 0.53, index = 2.79; chi2 = 1.32)  (small line)
6. (12) [okay?]     (nH = 0.59, index = 2.75; chi2 = 1.15)  (clean)

7. (18) [contam]    (nH = 0.52, index = 2.86; chi2 = 3.63)  (VISIBLE LINES)
8. (19) [contam]    (nH = 0.45, index = 2.54; chi2 = 2.49)  (VISIBLE LINES)
9. (21) [contam]    (nH = 0.68, index = 2.82; chi2 = 2.21)  (VISIBLE LINES)
10. (32) [mushy]    (nH = 0.9 , index = 2.98; chi2 = 1.13)  (clean)
11. (33) [mushy]    (nH = 0.8 , index = 2.91; chi2 = 1.15)  (clean, sparse)

12. (34) [okay?]    (nH = 0.68, index = 2.78; chi2 = 0.97)  (small line)
13. (35) [mushy]    (nH = 0.67, index = 3.00; chi2 = 3.40)  (VISIBLE LINES)
14. (36) [okay?]    (nH = 0.69, index = 2.85; chi2 = 1.11)  (clean)
15. (37) [okay?]    (nH = 0.76, index = 2.85; chi2 = 1.06)  (clean)
16. (38) [GOOD]     (nH = 0.72, index = 2.91; chi2 = 1.11)  (very small Si ln)
17. (40) [GOOD]     (nH = 0.72, index = 2.79; chi2 = 1.22)  (very small S line)

Results: fiddle with 1,5,12,16,17.  Discard 11.  Reject 7,8,9,13.
This leaves us: 1,2,3,4,5,6,10,12,14,15,16,17 -- 12 distinct regions, from
approx. 5 distinct filaments (counting the SW area into two).

Remark -- optimization of the regions (alignment / maximizing height /
minimizing width) is not too helpful, because we also have to address the
spectra as well.

All the spectra have a little bit of line present, most obvious in the
residuals (less so elsewhere). Both the Silicon XII line (~1.8 keV) and
Sulfur XV line (~2.4 keV) are usually, faintly visible.

Region 17(40) in particular has a more visible sulfur line.


Regions (good-3)
================

Discarded regions 7 (18), 8 (19), 9 (21), 13 (35)
Twiddled with the following regions.  Again looking at RED band profiles.

Region 1 (1): tried to cut back slightly farther? Narrowed slightly to avoid
              thermal knot to north
Region 5 (11): moved slightly down and rotated (avoid thermal towards north)
Region 12 (34): moved slightly down and rotated (avoid thermal towards north,
                where bad region 13/35 used to be)
Region 16 (38): narrowed to avoid thermal to north / narrow up the peak
                cut back a little more to avoid thermal stuff in back,
                slight rotation (less than 0.3 degrees)
Region 17 (40): narrowed region to avoid thermal to south, and moved more north
                slight rotation as well (to realign)

Saved new set of regions as `regions-good-3.reg`
Began the processing chain anew, regions look pretty good.


Regions (good-3-allback)
========================

Generating "addback" regions (vs. cutback...).  This is to give additional
downstream data for the profile fitting.  May introduce slight rotation to
improve fit quality =/...

Regions could stand to be twiddled a little more?
(numbers using regions-all numbering)
* Region 3 could be extended /drawn farther, and moved up a tiny bit
* Region 9 (4) could be split into two, narrowed up, twiddled with.
  It also benefits from a bit of rotation (try 47.15 degrees)
* Region 11 (5) is crap.  Twiddle that... but I don't know how to get around
  the faintness and multiple emission lines.  Maybe look at the spectra around
  there...
* Region 32 (7) is kinda iffy in red... hard to get a FWHM out.
  Same problem of mushiness, just a pile of mess.
  Speculatively, I wonder if messing with the energy bands would help...
* Region 37 (10) is also super ratty.  I have twiddled it alot in this addback
  iteration... copy this over to `regions-4-good` and twiddle more...

At least region 38 (11) is pretty nice...


Regions (good-4-ext) (2014 June 26-27)
======================================

## Priorities
1. Per Brian's suggestion, keep regions that look bad in 0.7-1 keV if they
   look okay in other bands (check 1-2 keV for shape, 4.5-7 keV for counts).
2. Don't extend regions too far ahead of shock; it will skew FWHM fits.
   Fewer data points will 1) improve fit of FWHM peak, 2) help shrink FWHM
   uncertainty (stretching of peak will affect chi-squared more)
3. Choose regions WITH thermal uptick in back.  Verify that they may be fit
   with FWHMs, then generate regions WITHOUT thermal uptick.  Then, further
   remove regions with contamination.
   (note -- to generate regions with/without upticks, do it by hand -- just
   try to keep the angle constant)

## Some ways to discriminate good/bad regions
With profile fits + eqwidth calculations, we can give more quantitative
discriminants for region selection:

1. Can the FWHM be resolved in 0.7-1 keV band?
   Looking at `good-3-allback` -- this is the ONLY band in which the FWHM is
   NOT consistently resolved.  All other bands have clear FWHMs.
2. How many counts are there in 0.7-1 keV, or 4.5-7 keV?
   Need more counts to decrease FWHM uncertainty, and hence better constrain
   energy dependence of filament widths
3. What's the estimated equivalent width of the Si line?
   Generally, 0.1 keV or less is good. Current set of regions is okay....

Protip: disable the y-axis gridlines, to better gauge the relative difference
between peak and trough

## Region "good-4-ext" selection

I review good-3-allback according to the above criteria (FWHM quality and
uncertainty, number of counts, Si line width; all using 5-band split), and
adjust the regions to create "good-4-ext" ("ext" = "allback" != "cutback").
No new regions added.


* Region 1 (1) FWHM good in all bands
               want more counts in 0.7-1 keV / 4.5-7 keV
               eqwidth = 0.12 keV
               ADJUSTMENT: lengthened towards north.  Narrowing didn't help the
               contamination issue anyways, it might be due to the thermal knot
               to the south.  Rotated slightly.
* Region 2 (3) FWHM NOT resolved in 0.7-1 keV, close to limit in 1-2 keV
               eqwidth = 0.10 keV
               ADJUSTMENT: narrowed, improving 1-2 keV peak (FWHM definitely
               measurable, still near limit) (0.7-1keV still bad).
               Rotated slightly and moved down a bit.
* Region 3 (5) FWHM is BARELY resolved in 0.7-1 keV.  Passable.
               eqwidth = 0.12 keV
               try to lengthen integration along rim (notes on good-3-allback)
               ADJUSTMENT: moved slightly south to get more counts
               (barely improve ability to resolve 0.7-1 keV rim)

* Region 4 (9) FWHM good in all bands
               consider splitting into two regions
               eqwidth = 0.05 keV
               ADJUSTMENT: split into two regions.  Narrower regions okay due
               to high count rate.  Also helps improve peak shape (since there
               are two interacting filaments, but we're neglecting the dim one
               -- it's ~1 order of magnitude less intense).  The slightly
               distorted front of the peak may mess with the fit and hence FWHM
               estimate, slightly.  Need to check this in the finer bands.

                (added region 9.5 here)

* Region 5 (11) FWHM barely resolved in 0.7-1 keV (okay in others)
                want more counts in 0.7-1 keV / 4.5-7 keV
                eqwidth = 0.2 keV (!!!!!)
                This region also has a noticeable sulfur line
                ADJUSTMENT: widened integration length, moved slightly south of
                the bright knot (maybe that's source of contamination).
                Not usable for 0.7-1 keV, period.
* Region 6 (12) FWHM good in all but 0.7-1 keV
                eqwidth = 0.05 keV
                ADJUSTMENT: widened and moved to get more counts.  Still pretty
                bad in 0.7-1 keV, but hopefully can get a FWHM.

* Region 7 (32) FWHM very bad in 0.7-1 keV.  Peak shapes are broad.
                needs more counts in 4.5-7 keV.
                eqwidth = 0.07 keV
                ADJUSTMENT: widened to get more counts.  Peak shape is bad, but
                the shape is consistently weird.  0.7-1 keV is very difficult
                to improve, here I GIVE UP and this will only be usable above
                1 keV energy.
* Region 8 (34) FWHM not resolved in 0.7-1, close to limit in 1-2 keV.
                want more counts in 4.5-7 keV.
                eqwidth = 0.10 keV
                ADJUSTMENT: widened to get more counts.  Still very hard to get
                FWHM in 0.7-1, 1-2 keV bands.  Shape is funny but at least
                consistently so.
* Region 9 (36) FWHM contaminated in 0.7-1 keV, but good otherwise
                eqwidth = 0.05 keV
                ADJUSTMENT: widened very little, to get more counts.  Yes,
                0.7-1 keV FWHM is hard to resolve still.
* Region 10 (37) FWHM barely resolved in 0.7-1 keV
                 want more counts in 4.5-7 keV, but workable
                 eqwidth = 0.06 keV
                 ADJUSTMENT: widened where possible on either side... though
                 quality of 0.7-1keV is of course still iffy.
* Region 11 (38) FWHM well resolved in all bands
                 eqwidth = 0.09 keV
                 ADJUSTMENT: widened slightly towards north and south, just to
                 get more counts...
* Region 12 (40) FWHM well resolved in all bands
                 eqwidth = 0.02 keV (very good)
                 ADJUSTMENT: widened slightly towards north/south (by ~1 pixel)
                 to get more counts.

Basically, 0.7-1 keV looks bad.  Several regions need more counts (mainly to
help 4.5-7 keV).

For all regions, I shortened the upstream (ahead of shock) length.


Next set of ALL regions
=======================

how to review all regions?...  e.g. region 21 is terrible in 0.7-1 keV, but
actually usable at 1-2, 2-7 keV.

   plan: merge regions-good-4 with discarded regions from regions-good-2, and
   remaining regions in regions-all.  Create regions-all-5.
   This is ad hoc as hell... I guess as long as I list the specific
   criteria/goals I aim to hit with the regions.
