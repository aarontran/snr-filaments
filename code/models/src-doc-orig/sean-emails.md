Emails and comments from Sean
=============================

Reproduced from my email inbox, for reference.  Many are also cc'ed to Brian,
Rob if information needs to be scrounged.

Tuesday, 2014 July 01
=====================

This email accompanies the code (`FullEfflength.f`, `fglists.dat`) and guide
(`FWHMGuide.pdf`) that Sean first sent.

> Attached is main program (FullEfflength.f), which prompts the user for an
> input of `B_0`, `eta_2` (the Diffusion Coefficient in terms of the Bohm value
> at 2keV, and mu.  It produces several files which contain radial profiles at
> different photon energies, the FWHM at each photon energy, and a diagnostic
> file of the electron distribution.  It requires fglists.dat to run. There are
> several modifications you can make in terms of output.
> 
> I hope you know Fortran 77!  Unfortunately that's what my advisor knew at the
> time and it's what everything is written in.  If you don't have a compiler
> already, gfortran seems to works fine.
> 
> I'm sure you will have lots of questions as the code was written originally
> in several bits and pieces and I wasn't expecting anyone else to have to use
> it, so as soon as you have them, ask away!

Thursday, 2014 July 03
======================

Email in response to some questions about code inputs (`m_E`, FWHM values).
N.B., I think references to equation (16) should be to equation (6) instead.

With this email, Sean also sent `Widthfun.py`, which fits an approximate
width-energy relation (equation (6)) to measured data.  As he notes, this model
incorporates less physics than the code implementing the full solution to the
continuous energy-loss transport equation (12).

(hence the two tables in the paper, one giving analytic and one giving
numerical solution best-fit values for eta (diffusion coefficient) and B0
(magnetic field))

## Estimating a single `m_E` for all energies

> > For Tycho, we have FWHMs for 4 to 5 energy bands, so we may estimate `m_E`
> > for adjacent bands (getting 3--4 values of `m_E`), or fit a power law and
> > get one `m_E` for the whole energy range.
>
> The reason you don't want to fit a power law of FWHM vs. energy for multiple
> bands is precisely because, as you say, `m_E` depends on energy.  What you
> would want want to do is fit the FWHM to Equation (16), though this is only
> an approximation without a lot of the more detailed physics.  In principle we
> want to extend the code I sent you to do the same thing, essentially, now
> with a more complicated and non-analytic function (treating the code as a
> function).  However, if this proves too difficult, fitting one value of `m_E`
> for the whole range when doing numerical fits seems like a good fallback. 
 
## Use of FWHM vs. `m_E` in code fitting
 
> > 1. Per your guide, to find `B_0/eta` using 3 FWHM values, you had wanted to
> >    find the best fit output FWHMs for the measured FWHMs.  In this case, do
> >    you need/want to calculate `m_E` from the code?
> >
> > (and, for SN 1006, in principle would it be the same if you had fitted to
> > the FWHMs at 1 keV and 2 keV, instead of `m_E` and FWHM at 2 keV?
> > excepting the fact that manual fitting would probably be much harder)
>
> Right.  The parameter `m_E` is basically an alternate way of describing the
> FWHM in two different energy bins.  You can either fit (nu1,FWHM1) and
> (nu2,FWHM2) or just (nu2,FWHM2) and (nu2,`m_E`).  Either one will give the
> same result.  

> > 2. Did you guys ever fit all 3 FWHM values to get `m_E`? (would that have
> >    been useful, anywhere?)  As opposed to the `m_E` calculation for
> >    adjacent bands.
>
> As described above, when dealing with 2+ energy bands, we fit function (16)
> to the data, for `D \propto E^\mu`, though I forgot to include that bit of
> code in my last email!  Attached is a Python routine that fits the shape of
> the (nu,FWHM) to function (16) in the paper.  It gives the results and the
> standard deviations and can easily be extended to an arbitrary number of
> bands.  Do you know any Python? It's much nicer to work with than Fortran!

## Energy dependence of `m_E`

> > 3. Does energy-dependence of `m_E` matter?  Fig. 3 seems to suggest that
> >    `m_E` decreases with increasing energy (per `$D ~ E^\mu$`) -- which
> >    would render less useful any "global" fit for `m_E`.
>
> Yes!  In SN1006, with only three bands, we weren't able to tell much about
> the dependence of `m_E` on E, but with 4-5 bands in Tycho this information
> could be very useful

Wednesday 2014 July 16
======================

I emailed a few questions:

* You are using Simpson's rule to integrate along lines-of-sight (yielding
  intensity profiles)?
* To confirm, the normalization of synchrotron emissivity doesn't matter?  Just
  curious as I saw a factor sqrt(nu * B) which looks like it doesn't affect the
  FWHM as well, maybe needed for magnetic damping / plots.  (rederiving, I got
  prefactor `$c_3 / (4 * \sqrt{c_1})$` )
* The guide you sent mentioned that there was some code for the magnetic
  damping model.  I'm wondering if that was included? as I couldn't see any
  code with a spatially varying B field.
* Do you think it's worth porting the code to Python?  I'm hope it would ease
  the addition of a fitting algorithm; the main issue I could see is slower run
  times (maybe okay with f2py, cython, or similar).  Just wondering if there
  are other issues you may foresee, given that this is your code!

Sean's replies:

1. Yep!
2. Also yes. I think I just stopped adding constants halfway through and never
went back to delete them. Though there is a factor that does matter which is
multiplied by one of the electron distributions to make sure the two versions
are consistently normalized, so if you remove any constants you will also have
to change that number appropriately.
3. Yeah, that's a second version of the code that I can send you tomorrow.
4. I do.  That is actually how I tried to add the fitting algorithm in the
first place (using f2py), but it didn't work immediately and I didn't have time
to debug at the time.  The porting worked, but unfortunately I'm not sure if I
still have the remnants of it...I will check. Also, I did this with an older
version of the code, so it might just be easier to redo it on your end if your
familiar with f2py.

(2014 July 17, Sean sends magnetic damping code `FullefflengthPohlab.f`)

Friday 2014 July 18
===================

## Sean reviewed some changes and comments to the code.

> To explain xmax and xmin:<br> xmax is the maximum distance away front the
> shock while xmin is the minimum distance away front the shock.  So xmax is
> actually to the left while xmin is to the right, to make things confusing.
> Also, the way the onedint array is indexed is so that i=1 corresponds to the
> shock and increasing i moves you farther away from the shock, and to top that
> off, delr(inu) is negative.

## One typo (or, set of typos) in electron distribution functions:

> > in the electron distribution functions for mu != 1  (equations 15, 16), I
> > am seeing:
> >
> >   `$\exp[ - n^{ 1 / (1 - \mu) } E / E_cut + ... ]$` in the paper, and
> >   `$\exp[ - n^{ -1 / (1 - \mu) } E / E_cut + ... ]$` in the code
> >
> > I couldn't find equations to cross-check in Lerche & Schlickeiser (1980),
> > since their electron spectrum didn't include an exponential cut-off (I
> > think).  Do you know which one is correct?
>
> Not another typo!  I believe the code version is correct.  The variable 'n'
> for mu>1 is basically the energy, and we want contributions from infinite
> energy to be zero, so the '-' sign ensures that the exponent (-1/(1-mu)) is
> positive and thus will rapidly increase as n goes to infinity so that the
> negative exponential will then be killed.  For mu<1 the variable n is
> basically 1/energy, so we want contributions to die off at n=0, which is
> ensured with the negative sign as the exponent (-1/(1-mu)) will be negative,
> so n raised to that power will blow up at zero and again kill the
> exponential. I will double check though.

Followup (dated Monday 2014 July 21):

> I double checked, and the code is indeed correct.  The variable change is `n
> = (E'/E)^(1-mu)`, with an integration variable energy `E'`, so that `E' = E
> n^(-1/(1-mu))` gets plugged into the source term of `exp(-E'/E_cut)`.  If I
> had to choose, I'm glad the code is the one that's correct!

## An explanation of the electron energy spectrum index s=2.2:

> The s=2.2 comes from the radio fit for the spectral index alpha = .6, which
> is defined such that the flux is propto \nu^{-alpha}.  You can relate alpha
> to s, where the electron energy spectrum is assumed to be \propto E^{-s}, by
> taking the synchrotron emissivity, making the integral dimensionless by
> converting the integration variable to x (defined after in equation 20), and
> then pulling out all the factors of /nu. You should get that the emissivity
> is proportional to nu^{(1-s)/2}.  They do it here
> http://www.cv.nrao.edu/course/astr534/SynchrotronSpectrum.html after equation
> 5C3 by taking the delta function approximation, but the approximation is not
> needed in order to get the scaling correct.

> s=2 is often taken as a general case as is derivable from the
> advection-diffusion in diffusive shock acceleration theory for some simple
> assumptions.  Other values might use different observational results or
> mightbe [sic] more sophisticated predictions obtained from diffusive shock
> acceleration simulations.


More?
=====

Further emails regarding the model/etc go here (possibly from Steve Reynolds as
well)


