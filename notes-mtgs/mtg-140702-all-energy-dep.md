Meeting to review some results
==============================

Looked over:
1. multiple FWHM calculations (as discussed with Brian)
2. global fit results for exponent `m_E` (using constant in power law vs. not?)
(holy shit it's 5th week)

To-dos
------
* Tabulate the results, like what Sean did in LaTeX -- put all the FWHM numbers
  together, with point-to-point `m_e` calculations and global fit `m_e`
* Send all the relevant files, plots, numbers, documentation.  Explain the
  procedure, how everything was done (for whatever paper in works).

* Profiles -- need to generate fits using intensity mosaics, because those have
  exposure and vignetting corrections.  Should make little difference though.
* Add more regions -- can we get the highest 2 energy bands, if the 2-3 keV
  band is going bad due to the sulfur line?
* Model code -- test it out first, before automating / wasting time scripting
  stuff up.  Parameter space could be weird, so might have trouble converging
* Why did they calculate it point to point in SN 1006?  Figure that out, maybe
  ask Sean.

Big question: what's the difference between calculating `m_E` pointwise, vs.
fitting to 3 points (or, 4-5 points for us)?  Why did Sean do this for SN 1006?
Email him and ask today (cc Rob and Brian)

If necessary to talk with Steve/Sean next week, avoid Tues/Weds (or make it
later) -- world cup semifinals for Rob.

* Kepler -- use CIAO `merge_obs` to get the mosaics Brian made, in flux units
  but don't get ahead, should finish the Tycho work first anyways.  But it does
  take a while to run.

X-ray imaging after Chandra?
----------------------------
A just curious question -- seems like it won't be for a while, given that
1. Chandra is still going pretty strong (even if its insulation is slowly
   degrading)
2. As Rob put it -- NASA went to the moon, so NASA can go to the moon, right?
   X-ray mirror manufacturers aren't around anymore?
   (also, obscenely expensive and obscenely heavy)
   Athena will be 5-10 arcsec... (follow-up to XMM-Newton)
