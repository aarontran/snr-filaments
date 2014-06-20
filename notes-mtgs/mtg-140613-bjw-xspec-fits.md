XSPEC fitting and stuff with Brian
==================================

Meeting ~2pm, Friday June 13, 2014

Start with the usual set-up.  We use the grouped spectrum, only occasionally
do we want to use ungrouped spectrum?  Remember to set `neivers`.

    cpd /xw
    data src2_grp.pi
    setplot energy
    ignore **-0.6,7.-**
    xset neivers 2.0

Now introduce our very, very simple model.  Simply multiply an absorption model
with a power law model.  Parameters are nH (column density, units 10^22 / cm^2),
PhoIndex (photon/spectral index), and norm (normalization).

    XSPEC12>model phabs * powerlaw
    ...
    1:phabs:nH>0.6
    2:powerlaw:PhoIndex>
    3:powerlaw:norm>
    ...

We look up an absorption for Tycho -- try nH = 0.6 (i.e., 0.6e22 cm^-2), from
fit values in Table 3 of Hwang, Decourchelle, Holt, Petre (2002, ApJ).
Let's freeze this in our fit.

    freeze 1

We can un-freeze it by:

    thaw 1
    fr 1  # re-freeze

To get information about our model, fit parameters...

    show model
    show parameters

Now let's go ahead and apply the fit.  Default is 10 iterations, if it hasn't
converged it will prompt you.

    fit 200 # number of iterations

Check that parameters are physically reasonable.  What's reasonable?
To give an idea, the data we plotted (very bright NW rim of Tycho), with the
absorption unfrozen, gave nH = 0.56, PhoIndex = 2.80, norm = 6.64e-4.
The reduced chi-squared was 1.18.  See some numbers in section below.

If the parameters look out of whack, we can explore the parameter space.
`steppar` will spew out a bunch of Chi-Squared values.
IF `steppar` finds a better parameter value (smaller chi-squared), it will ask
if you want to use that instead.

    steppar 2 1. 3. 50  # Vary PhoIndex from 1. to 3., with 50 steps

We can also set a parameter manually, although this will no longer be a "fit".

    newpar 1 0.6

Finally, we can plot our nice results using

    plot ldata residual

To interact with the plot (e.g., zoom in, save data to plaintext), invoke
interactive plotting mode with

    XSPEC12>ip
    PLT> rescale x 1.5 2.  # Zooming in on a faint Si line
    PLT> re y 0.01 .1
    PLT> exit
    XSPEC12>

Print out data.  There will be 4-5 columns: energy, energy error, flux,
flux error, and model values (but, verify this yourself).

    XSPEC12>iplot
    PLT> help WData
    PLT> wdata filename

Brian also noted that you can use `hardcopy` to print out plots, though I
mentioned here that I'd found `cpd file.ps/cps` to also work.  For reference,
usage of `hardcopy` appears to be:

    XSPEC12>hardcopy filename.ps color


XSPEC fits of filament spectra
------------------------------

For fits, email Brian/Rob if things come up.  If chi^2 is greater than ~1.5,
maybe run it by Rob or Brian.  Just run things by them if anything is odd.
Brian's suggested good numbers:

* Tycho absorption, expect nH ~ 0.4 to 0.8
* Power law index: ~2 to 3, 1 to 3 is also okay
* Normalization: whatever


Questions and feedback on backgrounds
-------------------------------------

### Q: Background selection -- does distance from SNR matter?  Because I see a
falloff in background counts farther away.

A: Yes, the chip counts degrade farther away from the center.  Could also be
the sheer brightness of Tycho.  But it's a small effect.

### Q: Background region shape -- can I work closer to the SNR, use a box?

A: It doesn't matter.  Can make an arbitrary polygon easily in DS9.

### Q: How do my backgrounds look? (shows Brian `profiles_all_bg.reg`)

A: Move the background regions a little farther out from the SNR to avoid
contamination / stray stuff.


Next steps + Nina's arrival
---------------------------

### Q: fitting radial profiles, are we just gonna use DS9's output?

A: Yes, it's super simple/straightforward.  Just averages pixels along
perpendicular lines.  Nothing to it.

### Q: fitting shape, do I have to write my own or can I just use SciPy etc?

A: Don't reinvent the wheel!  E.g., examples from the
[SciPy cookbook](http://wiki.scipy.org/Cookbook/FittingData).

Play with the models -- start from 1, 2 parameter models (e.g. a Gaussian),
then move forward.  See if simpler models work, whatever.

### Q: Next steps, after fits (hopefully next week)?

A: Remember we are investigating energy dependence of filament widths.
We could try a four band breakdown, since Tycho is so bright (2-4, 4-7 keV).
Of course this would use the same regions, just generate more plots/fits/data.

Look at the energy dependence.  In SN1006, the widths decreased with increasing
energy.  If we see the same thing here, that's interesting.  If not, that's also
interesting!  We call up Steve Reynolds, Sean Ressler, who agreed to help throw
this data into their model codes, see what comes out.

(Steve was Brian's Ph.D advisor; Sean coded up the models from Steve. Brian
asked if I'd be back at Berkeley, maybe could collaborate with Sean)


Nina's arrival next Monday
--------------------------

* Help get her set up with same software (ds9, CIAO, xtools)
* Download the same Tycho observations (Hughes' 750 ks)
* Start the `chandra_repro` etc. as she will need to process them all

