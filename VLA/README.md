README - Tycho VLA data
=======================

* `TYCHO_IF1.fits` -- VLA radio image, 1994, from Stephen Reynolds (in
  turn, from David Moffett?).  Presumably project AM0437.
  From header: units Jy/beam, freq 1.375 GHz.  Not sure if merged all of the
  separate VLA configurations' observations (all of A/B/C/D) or whatever.
  Also not sure if 1.375, 1.625 GHz freq observations could have been somehow
  combined or anything.

* `flmts.reg` (and `flmts.jpeg`) -- selection of 4 regions on Tycho's rims to
  look at varying X-ray and radio rim morphology, together

* `prf-flmts/` -- output from `code/regions/ds9projplotter.py`.  Profile data
  and plots from 1.375 GHz, 1-1.7 keV, and 4-7 keV

* `plot_radio_xray_prfs.py` -- one-off script to plots `flmts.reg` profiles
