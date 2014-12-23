README - TDEM0020, Pannuti Tycho observation
============================================

Here's the listobs output for reference.
For readability, I abbreviated some of the rows of text.

Questions for Jack
------------------

* When to worry about opacity / gain curve (if any considerations, other than
  just "at high freq"?
* How to inspect flags?  Check that online flagging hit all the shadowed
  antennas, etc?  Check which spectra/data/baselines/antenna are hit by flags?
  (partial answer -- table is too damn big to inspect by hand.  To guess which
  antennas would be shadowed, use some trigonometry...)
  (answer: `flagcmd`, see E. Momjian's slides from 2014 reduction workshop)

    flagcmd(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            inpmode='xml', action='plot', flagbackup=False)
  (`flagInfo = flagdata(vis='G55.7+3.4_10s.ms', mode='summary')` also works)
* RFI -- tfcrop vs rflag vs ??? (any general guidance from experience? I note
  they say rflag is better for broadband stuff, tfcrop for spiky stuff)
* What exactly is the purpose of this preliminary bandpass calibration?
  Why wouldn't we, e.g., first remove the most obvious RFI (using, say,
  tfcrop), then generate a bandpass to remove that effect, then run another RFI
  removal, then do the calibration?  I'm kind of confused as to what's going on
  / what's the intent of all this.
* When using `flagdata` -- it looks to me like it's already applying the same
  flags on all polarizations.  Do we still have to run the extension?
  (answer: in `tfcrop`, having extendflags=True does this already.  Same for
  `rflag`!  So no need to run again unless we want to extend even more
  aggressively, compared to defaults)

* How to identify and remove crosstalk?

* for myself- - what's the purpose of solving for antenna based delays again,
  and how exactly does this work?  conceptually, at least.

* spw:16,24 show strange increasing slope (odd bandpass?).  is this sufficient
  rationale to discard those windows?


* Does it make it harder for `tfcrop` to pick up the lines, or spread from the
  lines (e.g., if I flag 0:16 but some baselines require 0:15~17 to be
  flagged?)?  If the spectral windowing/whatever does consider the flagged data
  too, to cull things, then we're okay and flagging cannot do any harm.
  Otherwise, if order matters, we'd have to be careful...

* Preliminary calibration -- how do you know you aren't throwing RFI into the
  bandpass shape too?  How to best avoid that?!

Maybe useful: [AIPS walkthrough](http://www-astro.physics.ox.ac.uk/~hrk/AIPS_TUTORIAL/HRK_AIPS_1.html)
Just to see things to look out for!...

* General idea -- can we be more aggressive in applying algorithms to LR/RL
  correlation data, and then just extend them to the LL/RR correlation?

Some useful info
----------------

This is a 2 hour observation, looks to be divided into two 1 hr blocks for each
set of spectral windows.  Total 1 hr on target, 1 hr on calibrators, across 32
spectral bands (each with 64 channels?!)

16 spectral windows span 1008-1968 MHz with 64 MHz spacing, 64 ch (1 MHz) each
16 spectral windows span 1988-3884 MHz with 128 MHz spacing, 64 ch (2 MHz) each
So this solidly spans L and S bands.

In DnC configuration (fairly compact), expect HPBW between 7-46 arcsec.
Conveniently, largest angular scale resolvable is ~8 arcmin in S band so no
need to mosaic; remnant fits nicely in one pointing.
In S band near Tycho's declination (50-65 deg), RFI from Sirius XM may be
problematic?

Bandpass/flux calibrator: 0137+331 (3C48) (field 0)
Polarization angle/leakage calibrator: J0319+4130 (field 1)
Phase/amplitude calibrator: J2350+6440 (this is the important one) (field 2)
Tycho SNR (field 3)

Comparing to G55.7+3.4 tutorial, they used:
J1925+2106 is phase calibrator
3C147 is amplitude, bandpass calibrator



Reduction procedure?
--------------------

### plotants, observation log notes
Northern? arm extended (C-like?), SW/SE arms are more compact (D-like?)
ea12 is at edge of north; ea23 and ea10 are at edges of SW/SE respectively.
Central antennas are ea19, ea11, ea04, ea22, ea28.

Antennas 1, 2, 7, 8, 12 (these correspond to eaXX numbers, painted on dishes)
No other problems noted.  D config shadowing as usual (but I think the high
declination helps?).

Unless otherwise noted I'll use ea19 as the reference antenna.

### plotms notes

Plotting amplitude vs. freq for a few scans, I see strong RFI throughout!
[Catalog of sources](https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/rfi)

__L-band__ is a slew of mess.  I specify spwIDs in XX:XX, referring to
spectral window : channel (following CASA notation):

* 8:37 (1.557 GHz) - strongest peak, GPS
* 3:28, 3:46, 9:17 (1.228, 1.246, 1.601 GHz) - GPS, GLONASS (Russian GPS)
* 2:39~43 (~1.176 GHz) - broad peak, airplane navigation
* 14:27, 14:31, 15:17 (1.931, 1.935, 1.985 GHz) - cell phones

__S-band__ also strongly affected, but RFI more specific / localized, not as
many individual peaks/mess.

* 18:39 (2.322 GHz) - strongest, w/spillover/peaks to side.  Siriux/XM radio
* Generally, windows 17~19 are most contaminated

### Antenna position calibration

Run the command:

    >> gencal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
              caltable='TDEM0020.antpos', caltype='antpos')

     Determine antenna position offests from the baseline correction database
    offsets for antenna ea01 : -0.00660   0.00210   0.00100
    offsets for antenna ea02 : -0.00880   0.00200   0.00110
    offsets for antenna ea04 : -0.00150   0.00270   0.00150
    offsets for antenna ea06 : -0.00220   0.00240   0.00200
    offsets for antenna ea07 : -0.00750   0.00180   0.00180
    offsets for antenna ea11 : -0.00100   0.00000   0.00000
    offsets for antenna ea12 : -0.00940   0.00000   0.00000
    offsets for antenna ea13 :  0.00000   0.00070   0.00000
    offsets for antenna ea14 : -0.00170   0.00270   0.00000
    offsets for antenna ea15 : -0.00100   0.00450   0.00000
    offsets for antenna ea16 :  0.00150   0.00270   0.00000
    offsets for antenna ea17 : -0.00210   0.00240   0.00130
    offsets for antenna ea20 : -0.00210   0.00310   0.00190
    offsets for antenna ea22 : -0.00200   0.00250   0.00210
    offsets for antenna ea27 : -0.00260   0.00220   0.00240

     1 ~ 0.007 *
     2 ~ 0.009 *
     4 ~ 0.003
     6 ~ 0.002
     7 ~ 0.008 *
    11 ~ 0.001
    12 ~ 0.009 *
    13 ~ 0.001
    14 ~ 0.003
    15 ~ 0.005
    16 ~ 0.003
    17 ~ 0.002
    20 ~ 0.003
    22 ~ 0.003
    27 ~ 0.003

Where is antenna 8?  Looking at the EVLA baseline correction database
([link](http://www.vla.nrao.edu/astro/archive/baselines/)), it appears that
the corrections for antennas 1, 2, 7, 8, 12 were entered on Sep 21 -- before
this observation (ah, as noted in the log).  So those corrections are already
entered, but it just warns to keep an eye out.

Antennas 1, 2, 7, 12 were not moved again after Sept. 17, 2014.  So, on Oct. 12
they were able to determine a bunch more corrections (and for other antennas as
well).  These corrections are smaller than those just after the initial move
(~0.01 to 0.02 meter, vs. < 0.01 meter for the worst single-coord correction
here).  But, these are still substantial corrections.

Antenna 8 was moved very soon after our observation, so we don't have better
corrections for its position (as compared to other antennas).

Thus, we should keep an eye out for antenna 8's data quality!

### Atmospheric opacity, gain curve, stuff?

Working at low frequency so we are okay (hopefully).

### Flag antenna shadowing, zero-amplitude data

I think I requested online flagging to be applied, so the shadowing /
zero-amplitude data should be flagged.  Nevertheless, let's try running this
and see what happens.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='shadow',flagbackup=True)
    flagdata(vis='G55.7+3.4_10s.ms', mode='clip',
             clipzeros=True, flagbackup=True)

This takes a while (several minutes).

RFI removal
-----------

### Preliminary bandpass calibration and flagging

First apply Hanning smoothing (this modifies data column and is irreversible):

    hanningsmooth(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms', datacolumn='data')

Comparison images generated with:

    plotms(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
           scan='6', antenna='ea19', spw='', xaxis='freq', yaxis='amp',
           coloraxis='spw', symbolshape='circle', correlation='RR,LL',
           plotfile='amp_v_freq.beforeHanning.png') 
    plotms(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
           scan='6', antenna='ea19', spw='', xaxis='freq', yaxis='amp',
           coloraxis='spw', symbolshape='circle', correlation='RR,LL',
           plotfile='amp_v_freq.afterHanning.png')

    # Need to set y-axis to maximum at 20, else spw:18 dominates
plotms(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
       scan='19', antenna='ea19', spw='', xaxis='freq', yaxis='amp',
       coloraxis='spw', symbolshape='circle', correlation='RR,LL',
       plotfile='amp_v_freq.afterHanning_pt2.png',
       plotrange=[-1,-1,0,20])

Next, use phase calibrator to generate approximate bandpass calibration
Go through each spw and select RFI-free channels.  I copy the tutorial in
choosing 4 channels in each spectral window:

    plotms(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
           scan='6,10,17,21', antenna='ea19', xaxis='channel', yaxis='amp',
           iteraxis='spw', yselfscale=True, correlation='RR,LL', symbolshape='circle')

    0:30~33,
    1:30~33,
    2:17~20,
    3:5~8,      discard entirely (GPS/GLONASS) -- looks horrific
    4:38~41,
    5:30~33,
    6:30~33,
    7:30~33,
    8:8~11,     discard entirely (GPS) -- looks horrific, though less spazzy than 3
    9:48~51,    discard entirely (GPS/GLONASS)
    10:20~23,
    11:30~33,
    12:30~33,
    13:30~33,
    14:17~20,
    15:33~36,
    16:27~30,   strange increasing slope (amp vs channel)?!
    17:23~26,   discard entirely, likely
    18:17~20,   discard entirely? but some channels look okay
    19:42~45,   (see comment for spw 18)
    20:20~23,
    21:30~33,
    22:37~40,
    23:30~33,
    24:30~33,   shows strange increasing slope
    25:30~33,
    26:30~33,
    27:17~20,   bumpy? not sure if normal
    28:30~33,
    29:27~30,   weak RFI
    30:29~32,   weak RFI
    31:36~39    bumpy

Now write a gaincal command to get pre-bandpass phase calibration:

    spw_phase_cal=('0:30~33,1:30~33,2:17~20,3:5~8,4:38~41,5:30~33,'
        '6:30~33,7:30~33,8:8~11,9:48~51,10:20~23,11:30~33,12:30~33,13:30~33,'
        '14:17~20,15:33~36,16:27~30,17:23~26,18:17~20,19:42~45,20:20~23,'
        '21:30~33,22:37~40,23:30~33,24:30~33,25:30~33,26:30~33,27:17~20,'
        '28:30~33,29:27~30,30:29~32,31:36~39')

    gaincal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            caltable='TDEM0020.initPh',
            intent='CALIBRATE_PHASE*', solint='int', spw=spw_phase_cal,
            refant='ea19',
            minblperant=3, minsnr=3.0, calmode='p',
            gaintable='TDEM0020.antpos')

Logger notes that refant ea25 was used -- why?  Looking at the antennas over
time, it doesn't appear that ea19 (ID=18) ever dropped out.  Maybe at edges of
scans, ea19 came online just a tiny bit late or something, I don't know.
Here's the gaincal logger output:

    At 2014/09/24/11:02:02.5 (Spw=0, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=1, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=2, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/10:47:57.5 (Spw=3, Fld=2), using refant ea28 (id=26) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=3, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=4, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=5, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=6, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=7, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=8, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:42.5 (Spw=8, Fld=2), using refant ea28 (id=26) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=9, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=10, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=11, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=12, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=13, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=14, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/11:02:02.5 (Spw=15, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=16, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=17, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=18, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=19, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=20, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=21, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=22, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=23, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=24, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=25, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=26, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=27, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=28, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=29, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=30, Fld=2), using refant ea25 (id=24) (alternate)
    At 2014/09/24/12:01:52.5 (Spw=31, Fld=2), using refant ea25 (id=24) (alternate)

Looking at `plotants`, ea25 seems to be an acceptable choice for a reference
antenna.  Now I'm not sure whether it alternated between ea19 and ea25, which
would make things more painful.  I guess, we'll take a look at the phase
solutions and see.
Answer: it used ea19 in most cases where possible.  So I guess that's okay?

Looking through phase solutions: spws 3,8 are utter crap.
spws 30, 31 show some wiggling but look okay... (at least coherent)

Generate bandpass solution:

    bandpass(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms', caltable='TDEM0020.initBP',
             intent='CALIBRATE_PHASE*', solint='inf', combine='scan',
             refant='ea19', minblperant=3, minsnr=10.0,
             gaintable=['TDEM0020.antpos', 'TDEM0020.initPh'],
             interp=['','nearest'], solnorm=False)

Inspect the bandpass amplitude solutions to see what to flag by hand:

* 3:17~39   (mostly already flagged.  edges of window look okay)
* 5:52~53
* 8         (all of it)
* 9:1~3
* 14:26~27
* 17:40~41
* 18        (what happened here?!)

I'm not sure what to do about the wiggliness and smaller/more subtle RFI, or
continuum RFI, but the hope is that that will get flagged downstream.
Anyways, do all that and then apply the preliminary calibration:

    plotcal(caltable='TDEM0020.initBP', xaxis='chan', yaxis='amp',
            iteration='spw', antenna='ea01,ea12,ea19,ea25,ea23,ea10')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             spw='3:17~39,5:52~53,8,9:1~3,14:26~27,17:40~41,18')
    applycal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             gaintable=['TDEM0020.antpos', 'TDEM0020.initBP'], calwt=False)

Oof, I think I did something wrong.  Some of the RFI is gone, but there are
still some VERY obvious spikes left.

### Automatic flagging (using rflag)

Iterate through each spectral window -- generate a set of flags, twiddling
until the plots look good.  Then, extend flags along polarization (flags in one
polarization product should be extended to all of LL, RR, LR, RL), and also try
frequency/time.  Requires a lot of user interaction.

    # Extend along time/freq (inspecting display first)
    flagdata(vis='G55.7+3.4_10s.ms', mode='extend', 
             spw='0', growtime=50.0, growfreq=90.0,
             action='calculate', display='data',
             flagbackup=False)
    flagdata(vis='G55.7+3.4_10s.ms', mode='extend', 
             spw='0', growtime=50.0, growfreq=90.0,
             action='apply', display='')

Some notes:
* spw 0: rflag is having a lot of trouble (bandpass soln is bad, I think)
  I will try using tfcrop first, since there are some clear line features,
  particularly in circular polarization!
  Attempting to extend in time/freq (w/ growtime=50, growfreq=90) doesn't seem
  to do much.  I think `tfcrop` covers this pretty well already, in comparison
  to `rflag`.  So I don't extend in time/freq, just across polarization

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='tfcrop', spw='0', datacolumn='corrected',
             action='apply', display='', freqcutoff=2.5, timecutoff=3)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='extend', spw='0', extendpols=True,
             action='apply', display='')

* spw 1: very strong signal around channel 10, should have been manually
  excised.  What I'll do is set a high cutoff and cut this strong one, then
  come back and re-run flagdata to get the smaller stuff.
  Actually, I don't know if that made any difference (are statistics computed
  on all data, or just unflagged data, or does it not matter if we only care
  about local RMS. oh well).  I use rflag and go very aggressive here, which
  kills a lot of signal too.  Then extend along polarization, time, freq.
  We basically lose most of spw 1 in this way.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='rflag', spw='1', datacolumn='corrected',
             action='apply', display='', freqdevscale=5, timedevscale=5)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='rflag', spw='1', datacolumn='corrected',
             action='apply', display='', freqdevscale=2, timedevscale=2)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='extend', spw='1', extendpols=True,
             action='apply', display='')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='extend', spw='1', growtime=50.0, growfreq=90.0,
             action='apply', display='')

* spw 2: huge blobby signal right in the middle, clear in all correlations.
  Again running an aggressive removal.  But we're still missing a good deal,
  it's not even aggressive enough!

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='rflag', spw='2', datacolumn='corrected',
             action='apply', display='', freqdevscale=2.5, timedevscale=2.5)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='extend', spw='2', extendpols=True,
             action='apply', display='')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='extend', spw='2', growtime=50.0, growfreq=90.0,
             action='apply', display='')

* spw 3~31: because this isn't worth spending much time on, compared to working
  on actual measurements/data.  What I should do in the future is run `tfcrop`
  first, then generate the calibration and run `rflag`.

  This does artificially remove some "natural"-looking noise (~10%, looking at
  spw 26), but probably better to cull more real noise.
  In comparison, spw 17 is horrendous.
  Well, let's try this routine.  Then, inspect and cull by hand further.
  Started at 6:17pm.

    cmdlist = ["mode='tfcrop' freqcutoff=2.5 timecutoff=3.5",
               "mode='rflag' freqdevscale=2.5 timedevscale=2.5",
               "mode='extend' extendpols=True growtime=50.0 growfreq=90.0"]
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='list', inpfile=cmdlist, datacolumn='corrected',
             action='apply', display='')

  Shoot -- it looks like it might not be using the corrected data columns.
  This may not end well.  Oh well.

Data reduction redux (round 2)
------------------------------

Because the first attempt at preliminary bandpass calibration appeared so
unsuccessful, I'm going to try again from scratch.  I remove the calibration
and flags.  Note that the Hanning smoothing is _not_ removed from the `data`
column, and is irreversible (I think?) short of re-downloading the MS from the
JVLA archive.

### Initial flags

    clearcal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms')
    flagmanager(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
                mode='restore', versionname='Original')

    # These should have been applied online, but just to be safe
    # The logger suggests that zero-flagging does pick up some ~2% of data
    # (assuming it's not counting already-flagged data)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='shadow',flagbackup=True)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='clip', clipzeros=True, flagbackup=True)

    # Quack 1st sample of each scan (interval 5s from listobs)
    # Gets ~3% of data
    # (!) quack might be a bit too aggressive, especially for the short
    # phase calibrator scans -- we lose ~18.2% of the data!
    # I have no idea whether that's typical or not
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='quack', quackinterval=5.0, quackmode='beg',
             quackincrement=False, flagbackup=True)

I tried testing `tfcrop` first, but it's a little inconsistent.
I will try nailing the strongest RFI by hand first.
(I am not sure whether it's better to let `tfcrop` run first, then do this to
catch what `tfcrop` missed (on the rationale that ensuring `tfcrop` has extra
information, would work better...).  OR, if we run this first, then let
`tfcrop` and `rflag` pick up the leftovers...)

### Manual RFI flagging, first round (spikey strong stuff)

Use broadband amplitude plot to get a sense of typical (and atypical
amplitudes), quickly see major peaks and bad windows, etc.  Then iterate on
spectral windows, with plotms and flagdata, to flag indiv channels

    plotms(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
           scan='6,19', antenna='ea19', xaxis='channel', yaxis='amp',
           coloraxis='antenna2', iteraxis='spw',
           symbolshape='circle', correlation='RL,LR')

Windows:channels to flag, with some notes

    0:16;38~40
    1:17~20
    2:2~4;12~14;29~31;32~50
    3                 # signals overpoweringly strong in entire window
    4:0~14;29~32
    5:3~11            # two time-pulsed lines here
    8~9               # same problem as in 3
    10:32~48          # Swath of 3-4 lines, often pulsed
    14:24~60
    15:2~28

    16:19~23;48~50;56~63  # very clear, but also time-localized.  turn-on/off
                          # a line around 25~30 pops up in a minority of
                          # baselines (but when it appears, it's STRONG)
                          # also, 16 shows the continuous increase in amp
                          # w/ channel, like in 24. strange
    17                # 21~30 okay, rest too bad to salvage
    18                # dominated by satellite radio
    19:55~63
    20:28~34          # some baselines appear untouched...
                      # sporadic pulsed emission around ch 40,50 as well
    21:55~58          # occasional pulsed line, otherwise very clean
    22:26~28          # one strong line, mainly

    27:26~33;49~51    # Appear only on occasion, but flag them all
    28                # at least ~6 different lines, but inconsistent.
                      # either okay, or terrible. flag whole window
    29                # see 28
    30                # see 28.  e.g.,
                      # ea04 && ea19, ea05 && ea19, ea09 && ea20 are bad.
    31                # see 28

Clean windows, but run `tfcrop` anyways: 6,7,11,12,13, 23,24,25,26.
I know a few baselines in 7 have lines, for sure.

flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
         spw=('0:16;38~40,1:17~20,2:2~4;12~14;29~31;32~50,3,4:0~14;29~32,'
              '5:3~11,8~9,10:32~48,14:24~60,15:2~28,'
              '16:19~23;48~50;56~63,17,18,19:55~63,20:28~34,'
              '21:55~58,22:26~28,27:26~33;49~51,28~31'))

Wow, that flagged 50% of the data...
Note that spw:3,8,9,17,18,28,29,30,31 are entirely removed.




### Automatic excision

Now I try again using the `tfcrop` algorithm, to get the most egregious RFI
before generating a preliminary bandpass calibration.
Also fitting time variation to a polynomial, since we haven't calibrated phase
yet.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='tfcrop', datacolumn='corrected', timefit='poly',
             action='calculate', display='both', flagbackup=False)
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='tfcrop', datacolumn='corrected', timefit='poly',
             action='apply', display='')


### Preliminary bandpass calibration

For lack of better understanding, I use the same channels I picked before for
the phase calibration.  I hope it's okay.. we'll see.

    spw_phase_cal=('0:30~33,1:30~33,2:17~20,3:5~8,4:38~41,5:30~33,'
        '6:30~33,7:30~33,8:8~11,9:48~51,10:20~23,11:30~33,12:30~33,13:30~33,'
        '14:17~20,15:33~36,16:27~30,17:23~26,18:17~20,19:42~45,20:20~23,'
        '21:30~33,22:37~40,23:30~33,24:30~33,25:30~33,26:30~33,27:17~20,'
        '28:30~33,29:27~30,30:29~32,31:36~39')

    gaincal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            caltable='TDEM0020.initPh',
            intent='CALIBRATE_PHASE*', solint='int', spw=spw_phase_cal,
            refant='ea19',
            minblperant=3, minsnr=3.0, calmode='p',
            gaintable='TDEM0020.antpos')

    bandpass(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms', caltable='TDEM0020.initBP',
             intent='CALIBRATE_PHASE*', solint='inf', combine='scan',
             refant='ea19', minblperant=3, minsnr=10.0,
             gaintable=['TDEM0020.antpos', 'TDEM0020.initPh'],
             interp=['','nearest'], solnorm=False)

    applycal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             gaintable=['TDEM0020.antpos', 'TDEM0020.initBP'], calwt=False)

Looking at the bandpass solutions, we seem awfully short on points (when
plotting in the time domain), though the bandpass solution looks good.
Phase looks okay? but I can't evaluate the time stability at all.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='summary')

I see that about 42% of the data is flagged across all windows, fields, scans,
etc, on average.  (See `casapy-20141113-205245.log` for reference)

### Second round of automatic RFI removal...

This actually looks pretty helpful in extending flags along time, picking up
little extra tidbits, etc.  Also flag one more line by hand.  I'm too lazy to
inspect all the spws and scans manually, honestly.
Looking at a clean window (spw='11') it seems like we will lose some good data.
But, oh well.  We're just testing it out.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             spw='2:21~23')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             mode='rflag', datacolumn='corrected',
             timedevscale=2.5, freqdevscale=2.5,
             action='apply', display='')

Looking at just scans 6, 19; I'm seeing that 52% of the data is now flagged.
Okay.  Scan 1, spw='14:0~1;17~20;61' should be flagged, although scan 1 is so
short (and, being the first scan, maybe not too useful anyways).  In several
other scans I'm seeing stuff at spw='14:17~20' too.  S band looks fine.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             scan='1', spw='14:0~1;61')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             spw='14:17~20')

As suggested by the tutorial, I inspect some basic plots: amplitude vs.
frequency for several scans, amplitude vs. baseline for several spectral
windows, amplitude vs. UV distance.  UV dist plot I'm not quite sure what to
expect, though.  But anyways, all looks okay now.

Data calibration
----------------

(run 1x on 11/14, then run again on 11/16...)

### Flux/bandpass calibrator (3C48)

We pull up the model images in CASA for this.
L-band: spw='0~2,4~7,10~15' # 3,8,9 out
S-band: spw='16,19~27'      # 17,18,28,29,30,31 out

    setjy(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
          field='0137+331=3C48', model='3C48_L.im', spw='0~2,4~7,10~15',
          scalebychan=True, standard='Perley-Butler 2013')

    setjy(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
          field='0137+331=3C48', model='3C48_S.im', spw='16,19~27',
          scalebychan=True, standard='Perley-Butler 2013')

Now for phase calibration thingy.  This is from before (the first attempt at
RFI removal), but I remove the windows that have been flagged.  I pray that
these channels haven't been flagged much, but who knows...

Ah, because I flagged channels '14:17~20' just above, we need a new selection
here...  woops.  I use '14:10~13' now.

    spw_phase_cal=('0:30~33,1:30~33,2:17~20,4:38~41,5:30~33,6:30~33,7:30~33,'
                   '10:20~23,11:30~33,12:30~33,13:30~33,14:10~13,15:33~36,'
                   '16:27~30,19:42~45,20:20~23,21:30~33,22:37~40,23:30~33,'
                   '24:30~33,25:30~33,26:30~33,27:17~20')

    gaincal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            caltable='TDEM0020.initPh.2',
            intent='*CALIBRATE_BANDPASS*', solint='int', spw=spw_phase_cal,
            refant='ea19',
            minblperant=3, minsnr=3.0, calmode='p',
            gaintable='TDEM0020.antpos')

One warning here:

     Insufficient unflagged antennas to proceed with this solve.
       (time=2014/09/24/10:26:57.5 field=0 spw=15 chan=0)

Moving on, we now compute residual antenna-based delays.
Typical delays are -6 to 4 nanosec.

    gaincal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            caltable='TDEM0020.K0', intent='*CALIBRATE_BANDPASS*',
            spw='0~2,4~7,10~16,19~27', refant='ea19', solint='inf',
            gaintype='K', combine='scan', minsnr=3.0,
            gaintable=['TDEM0020.antpos','TDEM0020.initPh.2'])

To avoid the stupid warnings, explicitly ask for unflagged spectral windows
spw='0~2,4~7,10~16,19~27'

    bandpass(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             caltable='TDEM0020.BP', refant='ea19',
             intent='*CALIBRATE_BANDPASS*', spw='0~2,4~7,10~16,19~27',
             solint='inf', combine='scan', minblperant=3, minsnr=3.0,
             gaintable=['TDEM0020.antpos', 'TDEM0020.initPh.2', 'TDEM0020.K0'],
             interp=['','nearest', 'nearest'], solnorm=False)

Slew of warnings about insufficient unflagged antennas !!!... I'll ignore them
all.  Whee.  Quality of amp/phase bandpasses looks ehh, but I don't have any
idea of how to tell good/bad, really.

Arbitrarily, I will set minsnr=3.0 (instead of 10.0 (!) used in the G55.7+3.4
tutorial), hoping we'll get a better solution out without discarding so many
data.

Realization: all the flagged antennas now (for minsnr=3.0) are for field 0, at
very precise times -- 10:30:55.5 to 10:30:59.3, and 11:30:23.8 to 11:30:25.1.
Interesting.  Maybe we should just flag those entire time ranges.

    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             timerange='2014/09/24/10:30:55~2014/09/24/10:31:00')
    flagdata(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             timerange='2014/09/24/11:30:23~2014/09/24/11:30:26')

The attempt to flag the 2nd burst of bad data seems to be problematic.
I get a `MSSelectionNullSelection` exception, complaining that there are zero
rows...

### Per-antenna gain calibration

Use solint='inf' to get single solution for all scans.
Again, try with minsnr=3.0 instead of 10.0 to get more solutions?

    gaincal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
            caltable='TDEM0020.phaseAmp', intent='*PHASE*,*AMPLI*,*FLUX*',
            spw='0~2,4~7,10~16,19~27', solint='inf', refant='ea19',
            minblperant=3, minsnr=3.0,
            gaintable=['TDEM0020.antpos', 'TDEM0020.K0', 'TDEM0020.BP'])

Now get flux for the phase calibrator (J2350+6440, 3C468.1)

    fluxscale(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
              caltable='TDEM0020.phaseAmp',
              fluxtable='TDEM0020.phaseAmp.fScale',
              reference='0137+331=3C48', incremental=False)

Output values (range: 7.3 to 3.1 Jy over L band, decreasing) look reasonable
compared to VLA calibrator list data.  Correct magnitude and trend (decreasing
flux with increasing frequency).

### Apply calibration

Apply to calibrators for self-calibration, and to Tycho for science.
Note: I don't know why, but the G55.7+3.4 tutorial command does NOT include the
antenna delay calibration table (`*.K0`)

    applycal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             spw='0~2,4~7,10~16,19~27', intent='*TARGET*',
             gaintable=['TDEM0020.antpos', 'TDEM0020.K0',
                        'TDEM0020.BP', 'TDEM0020.phaseAmp.fScale'],
             calwt=False)
    applycal(vis='TDEM0020.sb29665410.eb29703711.56924.434018796295.ms',
             spw='0~2,4~7,10~16,19~27', intent='*PHASE*,*AMPLI*,*FLUX*',
             gaintable=['TDEM0020.antpos', 'TDEM0020.K0',
                        'TDEM0020.BP', 'TDEM0020.phaseAmp.fScale'],
             calwt=False, interp=['','nearest','nearest','nearest'])

Currently, the solutions look like crap.
* phases are centered on zero, but spread all the way around to +/- 180 deg.
  (solution: add 'TDEM0020.K0' in gaintable!)
* two balls offset slightly in amplitude, for phase calibrator.
  less evident when phases look good, BUT still clearly present.  It's not
  a visual artifact of an off scaling or anything.
* in both plots, there's a lot of mess near the bottom at amplitude 0-ish.

Field 2 (phase calibrator).
* spw:19, antenna ea15 is giving rise to the blobs off in phase by +/- 100 deg.
  Each blob corresponds to ea15 on one/other side of baseline, or something.
* spw:27, ea15 w/ the same blobs!  channel 4 gves rise to 8 distinct blobs.
* spw:0, channels 59, 22, 21 are giving rise to more blobs.
* Scans 4, 15 (the long scans on phase calibrator J2350+6440) are the gigantic
  blobs offset from the rest in amplitude, it seems.  And, these scans are
  associated with all the near-zero amplitude data, stuck at the bottom.

Flux calibrator has similar problems.
* spw:16,19,27 have weird stuff at 11:25:57.5.  Each blob is associated with a
  specific antenna (but not baseline).  Not surprising to have antenna
  association.  But the clustering in time seems to suggest RFI.
* When I plot data with no time/channel averaging, the bimodal amplitude blobs
  disappear for field 0, spw:1~2 only... maybe just a coincidence.

The big issue, then, is the amplitude calibration for everything.  Sigh.  Come
back to this...


### Multiscale, multi-frequency CLEAN


Listobs output
--------------

    ##########################################
    ##### Begin Task: listobs            #####
    listobs(vis="TDEM0020.sb29665410.eb29703711.56924.434018796295.ms",selectdata=True,spw="",field="",
            antenna="",uvrange="",timerange="",correlation="",scan="",
            intent="",feed="",array="",observation="",verbose=True,
            listfile="",listunfl=False,cachesize=50)
    ================================================================================
               MeasurementSet Name:  /Users/atran3/CASA-TDEM0020/TDEM0020.sb29665410.eb29703711.56924.434018796295.ms      MS Version 2
    ================================================================================
       Observer: Dr. Thomas G. Tom Pannuti PhD     Project: uid://evla/pdb/29665407  
    Observation: EVLA
    Data records: 7941024       Total integration time = 7165 seconds
       Observed from   24-Sep-2014/10:25:07.5   to   24-Sep-2014/12:24:32.5 (UTC)
       
       ObservationID = 0         ArrayID = 0
      Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent
      24-Sep-2014/10:25:05.0 - 10:25:55.0     1      0 0137+331=3C48            56160  [0, 1, ..., 15]  [5, 5, ..., 5] [UNSPECIFIED#UNSPECIFIED]
                  10:26:00.0 - 10:34:50.0     2      0 0137+331=3C48           595296  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_BANDPASS#UNSPECIFIED, CALIBRATE_FLUX#UNSPECIFIED]
                  10:34:55.0 - 10:42:50.0     3      1 J0319+4130              533520  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_POL_ANGLE#UNSPECIFIED, CALIBRATE_POL_LEAKAGE#UNSPECIFIED]
                  10:42:55.0 - 10:50:50.0     4      2 J2350+6440              533520  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  10:50:55.0 - 10:51:50.0     5      2 J2350+6440               61776  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  10:51:55.0 - 11:01:50.0     6      3 Tycho SNR               668304  [0, 1, ..., 15]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  11:01:55.0 - 11:02:50.0     7      2 J2350+6440               61776  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  11:02:55.0 - 11:12:45.0     8      3 Tycho SNR               662688  [0, 1, ..., 15]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  11:12:50.0 - 11:13:45.0     9      2 J2350+6440               61776  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  11:13:50.0 - 11:23:45.0    10      3 Tycho SNR               668304  [0, 1, ..., 15]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  11:23:50.0 - 11:24:45.0    11      2 J2350+6440               61776  [0, 1, ..., 15]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  11:24:50.0 - 11:25:45.0    12      0 0137+331=3C48            61776  [16, 17, ..., 31]  [5, 5, ..., 5] [UNSPECIFIED#UNSPECIFIED]
                  11:25:50.0 - 11:34:45.0    13      0 0137+331=3C48           600912  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_BANDPASS#UNSPECIFIED, CALIBRATE_FLUX#UNSPECIFIED]
                  11:34:45.0 - 11:42:40.0    14      1 J0319+4130              533520  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_POL_ANGLE#UNSPECIFIED, CALIBRATE_POL_LEAKAGE#UNSPECIFIED]
                  11:42:45.0 - 11:50:40.0    15      2 J2350+6440              533520  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  11:50:45.0 - 11:51:40.0    16      2 J2350+6440               61776  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  11:51:45.0 - 12:01:40.0    17      3 Tycho SNR               668304  [16, 17, ..., 31]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  12:01:45.0 - 12:02:40.0    18      2 J2350+6440               61776  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  12:02:45.0 - 12:12:35.0    19      3 Tycho SNR               662688  [16, 17, ..., 31]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  12:12:40.0 - 12:13:35.0    20      2 J2350+6440               61776  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
                  12:13:40.0 - 12:23:35.0    21      3 Tycho SNR               668304  [16, 17, ..., 31]  [5, 5, ..., 5] [OBSERVE_TARGET#UNSPECIFIED]
                  12:23:40.0 - 12:24:35.0    22      2 J2350+6440               61776  [16, 17, ..., 31]  [5, 5, ..., 5] [CALIBRATE_AMPLI#UNSPECIFIED, CALIBRATE_PHASE#UNSPECIFIED]
               (nRows = Total number of rows per scan) 
    Fields: 4
      ID   Code Name                RA               Decl           Epoch   SrcId      nRows
      0    NONE 0137+331=3C48       01:37:41.299431 +33.09.35.13299 J2000   0        1314144
      1    NONE J0319+4130          03:19:48.160102 +41.30.42.10305 J2000   1        1067040
      2    NONE J2350+6440          23:50:54.917601 +64.40.17.83000 J2000   2        1561248
      3    NONE Tycho SNR           00:25:14.000000 +64.08.39.00000 J2000   3        3998592
    Spectral Windows:  (32 unique spectral windows and 1 unique polarization setups)
      SpwID  Name           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) BBC Num  Corrs          
      0      EVLA_L#A0C0#0      64   TOPO    1008.000      1000.000     64000.0      12  RR  RL  LR  LL
      1      EVLA_L#A0C0#1      64   TOPO    1072.000      1000.000     64000.0      12  RR  RL  LR  LL
      2      EVLA_L#A0C0#2      64   TOPO    1136.000      1000.000     64000.0      12  RR  RL  LR  LL
      3      EVLA_L#A0C0#3      64   TOPO    1200.000      1000.000     64000.0      12  RR  RL  LR  LL
      4      EVLA_L#A0C0#4      64   TOPO    1264.000      1000.000     64000.0      12  RR  RL  LR  LL
      5      EVLA_L#A0C0#5      64   TOPO    1328.000      1000.000     64000.0      12  RR  RL  LR  LL
      6      EVLA_L#A0C0#6      64   TOPO    1392.000      1000.000     64000.0      12  RR  RL  LR  LL
      7      EVLA_L#A0C0#7      64   TOPO    1456.000      1000.000     64000.0      12  RR  RL  LR  LL
      8      EVLA_L#B0D0#8      64   TOPO    1520.000      1000.000     64000.0      15  RR  RL  LR  LL
      9      EVLA_L#B0D0#9      64   TOPO    1584.000      1000.000     64000.0      15  RR  RL  LR  LL
      10     EVLA_L#B0D0#10     64   TOPO    1648.000      1000.000     64000.0      15  RR  RL  LR  LL
      11     EVLA_L#B0D0#11     64   TOPO    1712.000      1000.000     64000.0      15  RR  RL  LR  LL
      12     EVLA_L#B0D0#12     64   TOPO    1776.000      1000.000     64000.0      15  RR  RL  LR  LL
      13     EVLA_L#B0D0#13     64   TOPO    1840.000      1000.000     64000.0      15  RR  RL  LR  LL
      14     EVLA_L#B0D0#14     64   TOPO    1904.000      1000.000     64000.0      15  RR  RL  LR  LL
      15     EVLA_L#B0D0#15     64   TOPO    1968.000      1000.000     64000.0      15  RR  RL  LR  LL
      16     EVLA_S#A0C0#16     64   TOPO    1988.000      2000.000    128000.0      12  RR  RL  LR  LL
      17     EVLA_S#A0C0#17     64   TOPO    2116.000      2000.000    128000.0      12  RR  RL  LR  LL
      18     EVLA_S#A0C0#18     64   TOPO    2244.000      2000.000    128000.0      12  RR  RL  LR  LL
      19     EVLA_S#A0C0#19     64   TOPO    2372.000      2000.000    128000.0      12  RR  RL  LR  LL
      20     EVLA_S#A0C0#20     64   TOPO    2500.000      2000.000    128000.0      12  RR  RL  LR  LL
      21     EVLA_S#A0C0#21     64   TOPO    2628.000      2000.000    128000.0      12  RR  RL  LR  LL
      22     EVLA_S#A0C0#22     64   TOPO    2756.000      2000.000    128000.0      12  RR  RL  LR  LL
      23     EVLA_S#A0C0#23     64   TOPO    2884.000      2000.000    128000.0      12  RR  RL  LR  LL
      24     EVLA_S#B0D0#24     64   TOPO    2988.000      2000.000    128000.0      15  RR  RL  LR  LL
      25     EVLA_S#B0D0#25     64   TOPO    3116.000      2000.000    128000.0      15  RR  RL  LR  LL
      26     EVLA_S#B0D0#26     64   TOPO    3244.000      2000.000    128000.0      15  RR  RL  LR  LL
      27     EVLA_S#B0D0#27     64   TOPO    3372.000      2000.000    128000.0      15  RR  RL  LR  LL
      28     EVLA_S#B0D0#28     64   TOPO    3500.000      2000.000    128000.0      15  RR  RL  LR  LL
      29     EVLA_S#B0D0#29     64   TOPO    3628.000      2000.000    128000.0      15  RR  RL  LR  LL
      30     EVLA_S#B0D0#30     64   TOPO    3756.000      2000.000    128000.0      15  RR  RL  LR  LL
      31     EVLA_S#B0D0#31     64   TOPO    3884.000      2000.000    128000.0      15  RR  RL  LR  LL
    Sources: 128
      ID   Name                SpwId RestFreq(MHz)  SysVel(km/s) 
      0    0137+331=3C48       0     -              -            
      0    0137+331=3C48       1     -              -            
      0    0137+331=3C48       2     -              -            
      0    0137+331=3C48       3     -              -            
      0    0137+331=3C48       4     -              -            
      0    0137+331=3C48       5     -              -            
      0    0137+331=3C48       6     -              -            
      0    0137+331=3C48       7     -              -            
      0    0137+331=3C48       8     -              -            
      0    0137+331=3C48       9     -              -            
      0    0137+331=3C48       10    -              -            
      0    0137+331=3C48       11    -              -            
      0    0137+331=3C48       12    -              -            
      0    0137+331=3C48       13    -              -            
      0    0137+331=3C48       14    -              -            
      0    0137+331=3C48       15    -              -            
      1    J0319+4130          0     -              -            
      1    J0319+4130          1     -              -            
      1    J0319+4130          2     -              -            
      1    J0319+4130          3     -              -            
      1    J0319+4130          4     -              -            
      1    J0319+4130          5     -              -            
      1    J0319+4130          6     -              -            
      1    J0319+4130          7     -              -            
      1    J0319+4130          8     -              -            
      1    J0319+4130          9     -              -            
      1    J0319+4130          10    -              -            
      1    J0319+4130          11    -              -            
      1    J0319+4130          12    -              -            
      1    J0319+4130          13    -              -            
      1    J0319+4130          14    -              -            
      1    J0319+4130          15    -              -            
      2    J2350+6440          0     -              -            
      2    J2350+6440          1     -              -            
      2    J2350+6440          2     -              -            
      2    J2350+6440          3     -              -            
      2    J2350+6440          4     -              -            
      2    J2350+6440          5     -              -            
      2    J2350+6440          6     -              -            
      2    J2350+6440          7     -              -            
      2    J2350+6440          8     -              -            
      2    J2350+6440          9     -              -            
      2    J2350+6440          10    -              -            
      2    J2350+6440          11    -              -            
      2    J2350+6440          12    -              -            
      2    J2350+6440          13    -              -            
      2    J2350+6440          14    -              -            
      2    J2350+6440          15    -              -            
      3    Tycho SNR           0     -              -            
      3    Tycho SNR           1     -              -            
      3    Tycho SNR           2     -              -            
      3    Tycho SNR           3     -              -            
      3    Tycho SNR           4     -              -            
      3    Tycho SNR           5     -              -            
      3    Tycho SNR           6     -              -            
      3    Tycho SNR           7     -              -            
      3    Tycho SNR           8     -              -            
      3    Tycho SNR           9     -              -            
      3    Tycho SNR           10    -              -            
      3    Tycho SNR           11    -              -            
      3    Tycho SNR           12    -              -            
      3    Tycho SNR           13    -              -            
      3    Tycho SNR           14    -              -            
      3    Tycho SNR           15    -              -            
      0    0137+331=3C48       16    -              -            
      0    0137+331=3C48       17    -              -            
      0    0137+331=3C48       18    -              -            
      0    0137+331=3C48       19    -              -            
      0    0137+331=3C48       20    -              -            
      0    0137+331=3C48       21    -              -            
      0    0137+331=3C48       22    -              -            
      0    0137+331=3C48       23    -              -            
      0    0137+331=3C48       24    -              -            
      0    0137+331=3C48       25    -              -            
      0    0137+331=3C48       26    -              -            
      0    0137+331=3C48       27    -              -            
      0    0137+331=3C48       28    -              -            
      0    0137+331=3C48       29    -              -            
      0    0137+331=3C48       30    -              -            
      0    0137+331=3C48       31    -              -            
      1    J0319+4130          16    -              -            
      1    J0319+4130          17    -              -            
      1    J0319+4130          18    -              -            
      1    J0319+4130          19    -              -            
      1    J0319+4130          20    -              -            
      1    J0319+4130          21    -              -            
      1    J0319+4130          22    -              -            
      1    J0319+4130          23    -              -            
      1    J0319+4130          24    -              -            
      1    J0319+4130          25    -              -            
      1    J0319+4130          26    -              -            
      1    J0319+4130          27    -              -            
      1    J0319+4130          28    -              -            
      1    J0319+4130          29    -              -            
      1    J0319+4130          30    -              -            
      1    J0319+4130          31    -              -            
      2    J2350+6440          16    -              -            
      2    J2350+6440          17    -              -            
      2    J2350+6440          18    -              -            
      2    J2350+6440          19    -              -            
      2    J2350+6440          20    -              -            
      2    J2350+6440          21    -              -            
      2    J2350+6440          22    -              -            
      2    J2350+6440          23    -              -            
      2    J2350+6440          24    -              -            
      2    J2350+6440          25    -              -            
      2    J2350+6440          26    -              -            
      2    J2350+6440          27    -              -            
      2    J2350+6440          28    -              -            
      2    J2350+6440          29    -              -            
      2    J2350+6440          30    -              -            
      2    J2350+6440          31    -              -            
      3    Tycho SNR           16    -              -            
      3    Tycho SNR           17    -              -            
      3    Tycho SNR           18    -              -            
      3    Tycho SNR           19    -              -            
      3    Tycho SNR           20    -              -            
      3    Tycho SNR           21    -              -            
      3    Tycho SNR           22    -              -            
      3    Tycho SNR           23    -              -            
      3    Tycho SNR           24    -              -            
      3    Tycho SNR           25    -              -            
      3    Tycho SNR           26    -              -            
      3    Tycho SNR           27    -              -            
      3    Tycho SNR           28    -              -            
      3    Tycho SNR           29    -              -            
      3    Tycho SNR           30    -              -            
      3    Tycho SNR           31    -              -            
    Antennas: 27:
      ID   Name  Station   Diam.    Long.         Lat.                Offset from array center (m)                ITRF Geocentric coordinates (m)        
                                                                         East         North     Elevation               x               y               z
      0    ea01  N12       25.0 m   -107.37.09.0  +33.54.30.0       -107.1669      870.2598       -7.3307 -1601110.034900 -5041488.083700  3555597.436200
      1    ea02  N16       25.0 m   -107.37.10.9  +33.54.48.0       -155.8413     1426.6388       -9.3855 -1601061.946800 -5041175.884000  3556058.032000
      2    ea03  E07       25.0 m   -107.36.52.4  +33.53.56.5        318.0421     -164.1874       -2.6819 -1600880.583700 -5042170.397400  3554741.463100
      3    ea04  E02       25.0 m   -107.37.04.4  +33.54.01.1          9.8166      -20.4427       -2.7746 -1601150.071100 -5042000.629400  3554860.721700
      4    ea05  W08       25.0 m   -107.37.21.6  +33.53.53.0       -432.1194     -272.1450       -1.5051 -1601614.093200 -5042001.650900  3554652.511800
      5    ea06  N06       25.0 m   -107.37.06.9  +33.54.10.3        -54.0677      263.8792       -4.2313 -1601162.592500 -5041828.994200  3555095.895300
      6    ea07  N14       25.0 m   -107.37.09.9  +33.54.38.5       -130.2486     1134.2144       -8.4821 -1601087.169300 -5041339.839500  3555815.854700
      7    ea08  N10       25.0 m   -107.37.08.2  +33.54.22.4        -86.6482      636.0548       -6.1235 -1601130.330000 -5041619.771300  3555403.733400
      8    ea09  W07       25.0 m   -107.37.18.4  +33.53.54.8       -349.9764     -216.7473       -1.7867 -1601526.378600 -5041996.850200  3554698.336400
      9    ea10  E09       25.0 m   -107.36.45.1  +33.53.53.6        506.0688     -251.8548       -3.5765 -1600715.938500 -5042273.189000  3554668.198000
      10   ea11  W02       25.0 m   -107.37.07.5  +33.54.00.9        -67.9879      -26.5288       -2.7337 -1601225.264000 -5041980.347590  3554855.693000
      11   ea12  N18       25.0 m   -107.37.12.0  +33.54.58.3       -183.8654     1746.5124      -10.3611 -1601034.381900 -5040996.527900  3556322.940300
      12   ea13  W04       25.0 m   -107.37.10.8  +33.53.59.1       -152.8647      -83.7939       -2.4621 -1601315.895500 -5041985.312070  3554808.313700
      13   ea14  E08       25.0 m   -107.36.48.9  +33.53.55.1        407.8311     -206.0306       -3.2217 -1600801.929000 -5042219.384400  3554706.431200
      14   ea15  E06       25.0 m   -107.36.55.6  +33.53.57.7        236.9059     -126.3503       -2.4648 -1600951.585000 -5042125.902000  3554772.989700
      15   ea16  W06       25.0 m   -107.37.15.6  +33.53.56.4       -275.8235     -166.7512       -2.0517 -1601447.195900 -5041992.513100  3554739.686600
      16   ea17  E04       25.0 m   -107.37.00.8  +33.53.59.7        102.7936      -63.7653       -2.6275 -1601068.804600 -5042051.916000  3554824.845500
      17   ea18  E03       25.0 m   -107.37.02.8  +33.54.00.5         50.6771      -39.4875       -2.7328 -1601114.352400 -5042023.153600  3554844.937600
      18   ea19  W01       25.0 m   -107.37.05.9  +33.54.00.5        -27.3531      -41.2994       -2.7502 -1601189.024540 -5042000.485700  3554843.424000
      19   ea20  N04       25.0 m   -107.37.06.5  +33.54.06.1        -42.6072      132.8482       -3.5378 -1601173.965600 -5041902.669500  3554987.527800
      20   ea21  E05       25.0 m   -107.36.58.4  +33.53.58.8        164.9878      -92.8081       -2.5370 -1601014.451700 -5042086.249300  3554800.790000
      21   ea22  N02       25.0 m   -107.37.06.2  +33.54.03.5        -35.6386       53.1794       -3.1447 -1601180.872480 -5041947.441700  3554921.622000
      22   ea23  W09       25.0 m   -107.37.25.2  +33.53.51.0       -521.9517     -332.7732       -1.2061 -1601710.024500 -5042006.915500  3554602.355100
      23   ea24  W05       25.0 m   -107.37.13.0  +33.53.57.8       -210.0976     -122.3871       -2.2585 -1601377.010700 -5041988.663500  3554776.394300
      24   ea25  W03       25.0 m   -107.37.08.9  +33.54.00.1       -105.3244      -51.7220       -2.5937 -1601265.137500 -5041982.549450  3554834.860400
      25   ea27  N08       25.0 m   -107.37.07.5  +33.54.15.8        -68.9061      433.1959       -5.0727 -1601147.936400 -5041733.823300  3555235.954800
      26   ea28  E01       25.0 m   -107.37.05.7  +33.53.59.2        -23.8729      -81.1415       -2.5834 -1601192.475300 -5042022.850400  3554810.447600
    ##### End Task: listobs              #####
    ##########################################

