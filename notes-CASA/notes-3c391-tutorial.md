README
======

Going through CASA tutorials for 3c391 and G55.7+3.4.

Idea: double check the high-level understanding with Jack...
(in particular, RFI removal for Tycho data since we are working at ~1 GHz


Questions
---------

* When are we supposed to use clearstat() exactly?
* In the 3C391 tutorial, what am I looking at in the final calibrated outputs?
How do I know this is good?  Especially for field 0 (where we're looking at the
cross-hand terms, LR/RL correlations).  For field 1 I can see that we have a
very stable calibration reference, at least.
* Corrected outputs -- what should I look for, esp. in the target data?
Tutorial suggests checking amplitude vs phase (noise should be random),
amplitude vs. uv distance (no RFI etc?)
* What are the main considerations when running clean (multiscale or w/e) -- it
  seems very subjective. When to stop? (check off-source RMS error,
  peak object brightness, ?)

> Note that the archive will automatically flag shadowed antennas as well as
> zero-valued data, if you request that online flags are applied.

* Why should zero-valued data be flagged?

* What's the difference between CLEAN cycles and iterations?  I think I saw
  this somewhere before..

* Converting Jy/beam to usable units to compare w/ X-ray data?

Useful links
------------

### Radio astronomy

    * [1](http://www.haystack.edu/edu/undergrad/materials/RA_tutorial.html)
    * [2](https://safe.nrao.edu/wiki/bin/view/Main/RadioTutorial)

### Calibration

    * [1](https://science.nrao.edu/science/meetings/2014/14th-synthesis-imaging-workshop/lectures-files/MoellenbrockCalibration2014_FINAL.pdf)
    * [2](https://science.nrao.edu/facilities/vla/ctw/drw2013/E_Momjian_Data_Reduction_Techniques_April2013.pdf)
    * [3](https://www.astron.nl/~mag/dokuwiki/lib/exe/fetch.php?media=radio_astronomy_lec_8_ma_garrett.pdf)

### Useful stuff

    * [1](https://science.nrao.edu/facilities/vla/docs/manuals/cal)

Calibration procedure summary
-----------------------------

### First review of data, antenna correction
* Inspect logs, listobs output, antenna configuration
* Select (good) antenna near center of array for calibration
* Flag data:
  - Bad antennas (log)
  - Dummy first scan
  - First sample from start of each scan
* Review data (`plotms`)
* Baseline correction, esp. for recently moved antennas (`gencal`)

### Visibility calibration
* Flux density (amplitude) calibration (`setjy`)
* Initial phase calibration (`gaincal`)
  - Plot phase over time for each antenna, ensure no jumps (`plotcal`)
* Delay calibration (`gaincal`)
  - Plot delays for each antenna (~ns) (`plotcal`)
* Bandpass calibration (`bandpass`)
  - Plot bandpass phase/amplitude over all channels and antennae (`plotcal`)
* Gain calibration (`gaincal`)
  - Previous calibrations were for flux density calibrator (J1331+3030)
  - Compute complex gains now for all calibrators
    (flux density, phase, polarization)
  - Plot gain phase/amplitude solutions
  - Check phase stability (difference of L/R polarizations)

### Polarization calibration
* Apply updated manual model w/ Stokes Q, U, V? (`setjy`)
  Check model spectrum (RR, RL amplitudes; RL phase)
* Cross-hand delays: `gaincal(gaintype='KCROSS', ...)`
* Leakage terms: `polcal(poltype='Df', ...)`.
  Check results with `plotcal`, expect leakages 5-15h
  Use `poltype='Df+QU'` if a calibrator polarization is unknown.  It solves for
  source polarization... then goes back to fill in leakage error???
* Position angle: `polcal(poltype='Xf', ...)`
  Check results with `plotcal`.  This is a correction for all antennae.
* Amplitude gains for secondary calibrators: `fluxscale`

### Apply calibration
* Apply tables to calibrators first (why? what does this do?) (`applycal`)
* Then, apply to science data (with linear interpolation) (`applycal`)
* Split off science data to new measurement set for convenience (`split`)

### Generate image (UV-plane to Stokes IQUV)
* Run multiscale CLEAN (`clean`), interactively


CASA tutorial on 3c391
----------------------

Supernova Remnant 3C391: 6cm Polarimetry and Continuum Imaging, Mosaicking

NOTE: I think my data may be corrupted, due to having ctrl-c'ed CLEAN several
times.

Workflow looks like:
* inspect logs, listobs output, antenna configuration
* flag any bad data, if known a priori (observing logs)
* `plotms` -- plot a slew of information
              E.g., visibility amplitude vs. baseline length
              (this is, roughly, a fourier transform of sky intensity)

* Baseline correction (recently moved antennas etc).
  Matches [database values](http://www.vla.nrao.edu/astro/archive/baselines/)
  I don't get how gencal is doing this; is it polling a table or computing
  corrections from the data?

    gencal(vis='3c391_ctm_mosaic_10s_spw0.ms',
           caltable='3c391_ctm_mosaic_10s_spw0.antpos',
           caltype='antpos')

* Obtain a flux density for the amplitude calibrator that was observed...
  "primary flux density calibrators" (3C138, 3C147, 3C286, 3C48 in the multiple
  EVLA bands), w/ "ultimate flux density scale" from WMAP observations of Mars?
  Ah -- WMAP goes down to K/Ka bands (20--40 GHz), which is covered by EVLA.

  In this example tutorial, our calibrators are:
  J1331+3030 = 3C 286 = visibility amplitude (flux) calibrator
    (also used for spectral bandpass, polarization position angle)
  J1822-0938 = visibility phase calibrator (closer in sky to target)
  J0319+4130 = 3C 84 = polarization calibrator
  (see equation relating observed and true visibilities)

    setjy(vis='3c391_ctm_mosaic_10s_spw0.ms', field='J1331+3030',
          standard='Perley-Butler 2010', model='3C286_C.im',
          usescratch=False,scalebychan=True,spw='')

* Initial phase calibration.  This will prevent decorrelation when data are
  vector averaged to get bandpass calibration solutions (?).

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.G0all', 
            field='0,1,9', refant='ea21', spw='0:27~36',
            gaintype='G',calmode='p', solint='int', 
            minsnr=5, gaintable=['3c391_ctm_mosaic_10s_spw0.antpos'])

  What's with the spectral windows? 0:27-36, how did they know to pick these?
  Okay, 0 is the spectral window (the only one). 27-36 are the channels.

  After removing antenna 5 (phase jumping around too much), we redo the phase
  calibration -- now without antenna 5, and for J1331+3030 only.

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.G0', 
            field='J1331+3030', refant='ea21', spw='0:27~36',
            calmode='p', solint='int', 
            minsnr=5, gaintable=['3c391_ctm_mosaic_10s_spw0.antpos'])

  What happened?  Earlier we solved for fields 0, 1, 9 (all 3 calibrators) and
  we got just the phase solution (calmode='p').  solint = 'int' indicates that
  we want solutions for 10 s intervals (the integration time) vs. over a longer
  time.  gaintype='G' is the default so not changed.

  Woops.  Remember to use tilde, not dash, for ranges.

    > *** Error *** MSSelection time error: Parse error at or near token '-'
    > (near char. 9 in string "08:02:00-08:17:00")
    > [TIP: Did you know we use "~" as the range operator (for a good reason)?]

* Delay calibration.  Gaintype 'K' gives delays.  Use solint='inf' to get one
  solution averaged over all times, scans (why?).  Use our antenna position +
  phase calibrations.

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.K0', 
            field='J1331+3030',refant='ea21',spw='0:5~58',
            gaintype='K', solint='inf', combine='scan', minsnr=5,
            gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                       '3c391_ctm_mosaic_10s_spw0.G0'])

  Resultant delays are a few nanosec.  refant `ea21` (index 18) of course has
  no delay in either L/R polarization

* Bandpass calibration.  Formerly skippable for continuum-only data.  But since
  EVLA data have finer spectral resolution (?), should apply to be safe, even
  if we only care about continuum data.

    bandpass(vis='3c391_ctm_mosaic_10s_spw0.ms',
             caltable='3c391_ctm_mosaic_10s_spw0.B0',
             field='J1331+3030',spw='',refant='ea21',solnorm=True,combine='scan', 
             solint='inf',bandtype='B',
             gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                        '3c391_ctm_mosaic_10s_spw0.G0',
                        '3c391_ctm_mosaic_10s_spw0.K0'])

  What is the "dynamic range" of spectral observations?
  Here we're using solint='inf', combine='scan' to merge all scans together,
  and compute correction for the whole observation.
  bandtype='B' indicates we solve channel by channel.  BPOLY allows polynomial
  fit to bandpass, experimental...

* Gain calibration (now, both amplitude and phase)
  First, compute solutions for the flux density calibrator

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.G1',
            field='J1331+3030',spw='0:5~58',
            solint='inf',refant='ea21',gaintype='G',calmode='ap',solnorm=F,
            gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                       '3c391_ctm_mosaic_10s_spw0.K0',
                       '3c391_ctm_mosaic_10s_spw0.B0'])

  spw='0:5~58'  again, ignore the edge channels affected by bandpass rolloff
  calmode='ap'  indicates both amplitude and phase
  solnorm=F     don't normalize solution amplitudes
  solint='inf'  one solution/scan (over all time)

  This one now replaces `3c391_ctm_mosaic_10s_spw0.G0`.
  Now, augment this table w/calibrations for phase and polarization calibrators

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.G1',
            field='J1822-0938',
            spw='0:5~58',solint='inf',refant='ea21',gaintype='G',calmode='ap',
            gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                       '3c391_ctm_mosaic_10s_spw0.K0',
                       '3c391_ctm_mosaic_10s_spw0.B0'],
            append=True)

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.G1',
            field='J0319+4130',
            spw='0:5~58',solint='inf',refant='ea21',gaintype='G',calmode='ap',
            gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                       '3c391_ctm_mosaic_10s_spw0.K0',
                       '3c391_ctm_mosaic_10s_spw0.B0'],
            append=True)

  Inspecting solutions, I see some drift in gain phases but it's smooth, so
  that's okay?  So what exactly is going on?

  > The complex gains for each antenna/spwid are determined from the data
  > column (raw data), divided by the model column, for the specified fields.
  > The gains can be obtained for a specified solution interval for each
  > spectral window, or by a spline fit to all spectral windows simultaneously. 

  Okay -- so we have NRAO's models for the calibrators, and we use them to
  scale our measurements (giving corrected/scaled measurements in time).
  How do we apply them to our measurements?

* Polarization calibration
  Going back to our initial run of `setjy` to calibrate amplitudes
  The stated intensities from `setjy` are:

  > J1331+3030 (fld ind 0) spw 0  [I=7.8169, Q=0, U=0, V=0] Jy,
  > (Perley-Butler 2010)

  > Scaling spw 0's model image by channel to I = [7.81908, 7.81688, 7.81468,
  > 7.81248, 7.81028, 7.80808, 7.80588, 7.80369, 7.80149, 7.7993, 7.79711,
  > 7.79492, 7.79273, 7.79055, 7.78836, 7.78618, 7.784, 7.78182, 7.77964,
  > 7.77746, 7.77528, 7.77311, 7.77093, 7.76876, 7.76659, 7.76442, 7.76225,
  > 7.76009, 7.75792, 7.75576, 7.7536, 7.75143, 7.74927, 7.74712, 7.74496,
  > 7.7428, 7.74065, 7.7385, 7.73635, 7.7342, 7.73205, 7.7299, 7.72776,
  > 7.72561, 7.72347, 7.72133, 7.71919, 7.71705, 7.71491, 7.71277, 7.71064,
  > 7.70851, 7.70637, 7.70424, 7.70211, 7.69999, 7.69786, 7.69574, 7.69361,
  > 7.69149, 7.68937, 7.68725, 7.68513, 7.68301] Jy (ch 0) for visibility
  > prediction.

  From I = 7.8169, fractional polarization 11.2% (C band around observation
  time), and polarization angle 66 degrees,
  compute the remaining Stokes parameters as:

    # math module included
    i0 = 7.81694 # Stokes I value for spw 0 ch 0
    p0 = 0.112*i0 # Fractional polarization=11.2%
    q0 = p0*cos(66*pi/180) # Stokes Q for spw 0
    u0 = p0*sin(66*pi/180) # Stokes U for spw 0

  We can also get a spectral index, using the two values farthest apart.  I am
  super confused as to the numbers they're using (I only see 64 channels with
  the last value being I=7.68301; frequencies 4536 MHz to 4664 MHz look right
  though).  I'll use the numbers I see in `setjy` to be sure:

    >>> alpha = log(7.68301/7.8169) / log(4664./4536.)
    -0.6208398712752974

  Eh, that's close enough, right?  Okay.
  Now we have to modify the measurement set headers, because `setjy` runs
  calibration for Stokes I only?  (why doesn't it just include Q, U, V too?)

    delmod('3c391_ctm_mosaic_10s_spw0.ms')
    # deletes "model visibility data representations"

    # PREVIOUSLY we used:
    # setjy(vis='3c391_ctm_mosaic_10s_spw0.ms', field='J1331+3030',
    #       standard='Perley-Butler 2010',
    #       spw='', model='3C286_C.im',
    #       scalebychan=True, usescratch=False)

    # THIS TIME we use:
    setjy(vis='3c391_ctm_mosaic_10s_spw0.ms', field='J1331+3030',
          standard='manual',
          spw='0', fluxdensity=[i0,q0,u0,0],
          spix=alpha, reffreq='4536.0MHz',
          scalebychan=True, usescratch=False)

  So we're supplying our own model now, with Stokes I, Q, U, and spectral index
  alpha.  The spectral index and I are derived from the NRAO model, anyways.
  But we now add our manually computed Q, U.

  Remember usescratch=False; do NOT write model visibilities to `.ms` file
  Note that we only have one spectral window (spw='' same as spw='0').
  If we have multiple windows, may need multiple models (manual/NRAO/whatever)

  __Solving for cross-hand delays__

    gaincal(vis='3c391_ctm_mosaic_10s_spw0.ms',
            caltable='3c391_ctm_mosaic_10s_spw0.Kcross',
            field='J1331+3030', spw='0:5~58',
            gaintype='KCROSS', solint='inf', combine='scan', refant='ea21',
            gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                       '3c391_ctm_mosaic_10s_spw0.K0',
                       '3c391_ctm_mosaic_10s_spw0.B0',
                       '3c391_ctm_mosaic_10s_spw0.G1'],
            gainfield=['','','','J1331+3030'],
            interp=['linear','nearest','nearest','linear'], parang=T)

  What's going on here?!
  Now that we have applied a polarization calibration model, we correct for:
  antenna position correction, RR/LL delays, bandpass, complex gain;
  then we compute a calibration for R vs. L delay

  > Note that if we did not solve for this delay, it would be absorbed into the
  > phases per channel of the following Df and Xf solutions. This would not
  > cause us problems, as we are not solving for the Q+iU polarization of our
  > D-term calibrator (we are using unpolarized 3C84 for that) but if we were
  > (e.g. using our gain calibrator J1822-0938 with poltype='Df+QU') then this
  > step would be essential. 

  _Aaron blinks_.  What?  Df/Xf solutions -- leakage and position angle
  calibrations for polarization, as a function of frequency (across channels)

  __Solving for leakage terms__

  Call with `poltype='Df'` to get leakage (D) on per channel (freq f) basis.
  Use unpolarized calibrator J0319+4130 (3C 84) to solve for the
  _instrumental polarization_.

    polcal(vis='3c391_ctm_mosaic_10s_spw0.ms',
           caltable='3c391_ctm_mosaic_10s_spw0.D1',
           field='J0319+4130',spw='0:5~58',
           refant='ea21',poltype='Df',solint='inf',combine='scan',
           gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                      '3c391_ctm_mosaic_10s_spw0.K0',
                      '3c391_ctm_mosaic_10s_spw0.B0',
                      '3c391_ctm_mosaic_10s_spw0.G1',
                      '3c391_ctm_mosaic_10s_spw0.Kcross'],
           gainfield=['','','','J0319+4130',''],
           interp=['linear','nearest','nearest','linear','nearest'])

  Okay.  I think I see what's going on...
  Linear interp for antenna position correction
  Nearest neighbor for delay correction
  Nearest neighbor for bandpass correction
  Linear interp for complex gain (computed on all our calibrators, but we
  select only the calibrator we're using)
  Nearest neighbor for polarization cross-hand (RL/LR) delays

  And we combine all scans together, integrate over all time, to get a single
  frequency dependent correction.

  Antenna ea04 (index 3) is missing.  Run polcal with `poltype='Df+QU'` to get
  a leakage solution for ea04.  Its parallactic angle coverage is good (?!?!)
  Solve simultaneously for source polarization.

  __Solving for polarization angle__

    polcal(vis='3c391_ctm_mosaic_10s_spw0.ms',
           caltable='3c391_ctm_mosaic_10s_spw0.X1',
           field='J1331+3030',combine='scan',
           poltype='Xf',solint='inf',
           gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                      '3c391_ctm_mosaic_10s_spw0.K0',
                      '3c391_ctm_mosaic_10s_spw0.B0',
                      '3c391_ctm_mosaic_10s_spw0.G1',
                      '3c391_ctm_mosaic_10s_spw0.Kcross',
                      '3c391_ctm_mosaic_10s_spw0.D2'],
           gainfield=['','','','J1331+3030','',''],
           interp=['linear','nearest','nearest','linear','nearest','nearest'])

  Note we are using a different source now!  Flux density calibrator gives
  position angle too

  __Secondary calibrators' flux__

  Unlike the primary flux calibrator, these fluxes aren't as well known and are
  time-dependent!

    Flux density for J1822-0938 in SpW=0 (freq=4.536e+09 Hz) is:
        2.34087 +/- 0.0063875 (SNR = 366.477, N = 46)
    Flux density for J0319+4130 in SpW=0 (freq=4.536e+09 Hz) is:
        13.9337 +/- 0.0386588 (SNR = 360.428, N = 44)

  Check the VLA calibrator
  [manual](https://science.nrao.edu/facilities/vla/docs/manuals/cal)
  to ensure fluxes look reasonable.  They vary (often being AGN...).

* Apply calibration
  First, apply to all calibrators (substitute names appropriately)

    applycal(vis='3c391_ctm_mosaic_10s_spw0.ms',
             field='J1331+3030',
             gaintable=['3c391_ctm_mosaic_10s_spw0.antpos', 
                        '3c391_ctm_mosaic_10s_spw0.fluxscale1',
                        '3c391_ctm_mosaic_10s_spw0.K0',
                        '3c391_ctm_mosaic_10s_spw0.B0',
                        '3c391_ctm_mosaic_10s_spw0.Kcross', 
                        '3c391_ctm_mosaic_10s_spw0.D2',
                        '3c391_ctm_mosaic_10s_spw0.X1'],
             gainfield=['','J1331+3030','','','','',''], 
             interp=['','nearest','','','','',''],
             calwt=[False],
             parang=True)

  We are inputting: antenna position, fluxscale (final version of complex gain
  tables: G0all, G0, G1), delay (RR/LL), bandpass, cross-hand delay (RL/LR),
  leakage correction (Df; here Df+QU to get antenna ea04), and position angle
  correction (Xf).  Okay, cool.

  calwt=[False], leave as is.  Doesn't work...
  parang=True, account for our polarization calibrations!

  NOW, finally, apply calibration to our actual target fields!
  We use the secondary calibrator J1822-0938; this is our visibility phase
  calibrator, which has been calibrated against the better characterized
  J1331+3030 for flux.  J1822-0938 is the calibrator closest in sky to our
  source.

    applycal(vis='3c391_ctm_mosaic_10s_spw0.ms',
             field='2~8',
             gaintable=['3c391_ctm_mosaic_10s_spw0.antpos', 
                        '3c391_ctm_mosaic_10s_spw0.fluxscale1',
                        '3c391_ctm_mosaic_10s_spw0.K0',
                        '3c391_ctm_mosaic_10s_spw0.B0',
                        '3c391_ctm_mosaic_10s_spw0.Kcross', 
                        '3c391_ctm_mosaic_10s_spw0.D2',
                        '3c391_ctm_mosaic_10s_spw0.X1'],
             gainfield=['','J1822-0938','','','','',''], 
             interp=['','linear','','','','',''],
             calwt=[False],
             parang=True)

  Now we use linear interp, to interpolate the calibration for our source
  observations.  Above, we want to use the nearest calibration (so that
  observations match up (?)).

* Split off science measurement set for further analysis

    split(vis='3c391_ctm_mosaic_10s_spw0.ms',
          outputvis='3c391_ctm_mosaic_spw0.ms',
          datacolumn='corrected',field='2~8')

* Imaging (CLEANing)
  First inspect amplitude vs. UV-distance.  Max baseline = 16000 wavelengths,
  angular scale is 12 arcseconds (1/16000) (recall synthesized beam is set by
  largest baseline).

  Best CLEANing with 3-5 pixels across synthesized beam (why?)
  Choose pixel size 2.5 arcsec
  Diameter is ~9 arcminutes = 540 arcsec = 216 pixels.
  Choose image size 480 x 480 (prefer composite number divisible by any pair of
  2,3,5 for FFTW).

    clean(vis='3c391_ctm_mosaic_spw0.ms',
        imagename='3c391_ctm_spw0_noms_I',
        field='',spw='',
        mode='mfs',
        niter=5000,
        gain=0.1, threshold='1.0mJy',
        psfmode='clark',
        imagermode='mosaic', ftmachine='mosaic',
        multiscale=[0], 
        interactive=True,
        imsize=[480,480], cell=['2.5arcsec','2.5arcsec'],
        stokes='I',
        weighting='briggs',robust=0.5,
        usescratch=False)

  Let's see...
  mode='mfs' indicates multi-frequency synthesis imaging.  account for
  corrections to uv-distance due to the finite bandwidth of our images
  (different frequencies -> different distances)
  niter/threshold determine when to stop CLEAN process
  interactive -- pick out polygons in viewer, and change cleaning process as
  you go (as you see where should not be cleaned when low-intensity emission
  pops out)
  briggs weighting -- correction for lack of short spacings (less sensitive to
  large scale features, more patchiness/noise)
  imagermode, ftmachine -- mosaic uses all fields' data.

  Seems to be hanging (interactive window doesn't appear, doesn't seem to do
  much?).  I killed CASA a few times, deleted and regenerated the split-off
  science data, and restarted the cleaning process.

  Answer: it just takes a while, deconvolution is hard.  Letting it run...

### Tutorial part 2 (multiscale clean, image analysis, self-calibration)

* Multiscale cleaning.  Run with the command:

    clean(vis='3c391_ctm_mosaic_spw0.ms',imagename='3c391_ctm_spw0_IQUV',
          field='',spw='',
          mode='mfs',
          niter=25000,
          gain=0.1, threshold='1.0mJy',
          psfmode='clarkstokes',
          imagermode='mosaic', ftmachine='mosaic',
          multiscale=[0, 6, 18, 54], smallscalebias=0.9,
          interactive=True,
          imsize=[480,480], cell=['2.5arcsec','2.5arcsec'],
          stokes='IQUV',
          weighting='briggs',robust=0.5,
          pbcor=False,
          usescratch=False)

  Now using 'clarkstokes' psfmode w/ 4 polarizations.  Multiscale clean uses
  components w/ sizes (pixels) 0, 6, 18, 54.  6 pixels is comparable to beam
  size. 
  pbcor False doesn't correct for primary beam (do later).
  But this means the output image is not flux correct.
  Run `immath` to divide out the primary beam response manually

    immath(outfile='3c391_ctm_spw0_IQUV.pbcorimage',
           mode='evalexpr',
           imagename=['3c391_ctm_spw0_IQUV.image','3c391_ctm_spw0_IQUV.flux'],
           expr="IM0[IM1>0.2]/IM1")

  NOTE: select "All polarizations" button BEFORE drawing mask in CASA viewer

  Use viewer to inspect images and compute on/off source statistics (mean flux,
  peak flux, RMS noise etc in Jy/beam).

* Construct polarization images from datacube (IQUV).  We can compute the total
  linear polarization (P = sqrt(Q^2 + U^2)) and the position angle, and get out
  polarization vectors.  Check V image in particular, circular polarization can
  indicate unflagged RFI (since it's not typically produced by astrophysical
  sources).

  Extract from flat noise images specifically:

    immath(imagename='3c391_ctm_spw0_IQUV.image',outfile='3c391_ctm_spw0.I',expr='IM0',stokes='I')
    immath(imagename='3c391_ctm_spw0_IQUV.image',outfile='3c391_ctm_spw0.Q',expr='IM0',stokes='Q')
    immath(imagename='3c391_ctm_spw0_IQUV.image',outfile='3c391_ctm_spw0.U',expr='IM0',stokes='U')
    immath(imagename='3c391_ctm_spw0_IQUV.image',outfile='3c391_ctm_spw0.V',expr='IM0',stokes='V')

  Form a linear polarization image (P^2 = sqrt(Q^2 + U^2)):

    immath(outfile='3c391_ctm_spw0.P_unbias',
           mode='poli',
           imagename=['3c391_ctm_spw0.Q','3c391_ctm_spw0.U'],
           sigma='0.000041Jy/beam')

  I don't see where they got sigma = 0.041 mJy/beam from?!
  sigma must be given in Jy/beam.
  Anyways, filter out NaNs (meant for CASA 4.1 or older, but to be safe...) and
  apply primary beam correction:

    immath(outfile='3c391_ctm_spw0.P_unbias_filtered',
           mode='evalexpr',
           imagename=['3c391_ctm_spw0.P_unbias'],
           expr="IM0[IM0>-10000.0]")

    immath(imagename='3c391_ctm_spw0_IQUV.flux',
           outfile='3c391_ctm_spw0.Qflux',
           expr='IM0',stokes='Q')
    immath(outfile='3c391_ctm_spw0.pbcorP',
           mode='evalexpr',
           imagename=['3c391_ctm_spw0.P_unbias_filtered','3c391_ctm_spw0.Qflux'],
           expr="IM0[IM1>0.2]/IM1")

  Here we extracted one plane of the primary beam flux file (same for all),
  then used the same corrrection as earlier to subtract out beam response.

  Generate a polarization angle image from Stokes Q, U.
  We set a threshold to avoid computing the angle for noise areas
  (recall off-source RMS was ~0.055 mJy/beam in Stokes Q,U images, no pbcor)

    immath(outfile='3c391_ctm_spw0.X',
           mode='pola',
           imagename=['3c391_ctm_spw0.Q','3c391_ctm_spw0.U'],
           polithresh='0.2mJy/beam')

  For bonus, generate a fractional linear polarization image too

    immath(outfile='3c391_ctm_spw0.F',
           mode='evalexpr',
           imagename=['3c391_ctm_spw0.I','3c391_ctm_spw0.Q','3c391_ctm_spw0.U'],
           expr='sqrt((IM1^2+IM2^2)/IM0[IM0>1.2e-3]^2)')

* Generate combined image -- load all three images
  (use LEL expression to only plot polarization vectors where linear
  polarization is sufficiently strong).
  Twiddle plot settings appropriately

* Spectral index imaging (exercise to the reader)

  Supply two images at different spectral frequencies and use `immath` to do
  the calculation.  This whole exercise was completed in one spectral window
  (spw0) as provided for download.  We could work using the highest and lowest
  channels -- how would this be done?

  The clean image we generate is synthesized at 4.6 GHz with a 128 MHz
  bandwidth (each channel is 2 MHz wide).  We used multi-frequency synthesis in
  clean to collapse all the channels together, accounting for the slight
  frequency differences in the UV plane to get a better result.

  Solution -- from the calibrated dataset, specify the spectral window +
  channels to use.  For us we'd want `spw='0:0'` and `spw='0:63'`.
  Let's give it a shot real quickly (I expect clean should run faster this
  way...)

    clean(vis='3c391_ctm_mosaic_spw0.ms',imagename='3c391_ctm_spw0_I_ch63',
          field='',spw='0:63',
          mode='mfs',
          niter=25000,
          gain=0.1, threshold='1.0mJy',
          psfmode='clark',
          imagermode='mosaic', ftmachine='mosaic',
          multiscale=[0, 6, 18, 54], smallscalebias=0.9,
          interactive=True,
          imsize=[480,480], cell=['2.5arcsec','2.5arcsec'],
          stokes='I',
          weighting='briggs',robust=0.5,
          pbcor=False,
          usescratch=False)

  So I run a multiscale clean on Stokes I only, for channels 0, 63.  Only run
  two polygon cleaning rounds on each channel (1000 iterations) to save time
  (so data quality will be bad but that's okay).

* Self-calibration

Here's the [tutorial section](http://casaguides.nrao.edu/index.php?title=EVLA_Advanced_Topics_3C391#Self-Calibration).
I'm passing over running this for now (especially given how slow cleaning is on
this machine), but conceptually it seems straight forward.  Take a look at the
VLA data reduction workshop materials.




