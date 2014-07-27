README
======
2014 July 26

Before reading code, please read `rminarc-bug.md` and/or talk to Aaron about
things that may come up.  The main script you need is `models_all.py`, and
tables which may have been already generated in `tables/`.


Main model code
---------------
Modified version of Sean's code.  See writeup on Sean's code, which checks and
explains details and corrections to Sean's original code.
The backup version dated 2014 July 21 resembles Sean's code in its output and
layout.  The newer version is more aggressively edited to use with f2py.

    old-code/FullEfflength_mod_2014-07-21_working.f
    FullEfflength_mod.f
    fullmodel_recompile.py  # Code snippet to recompile FullEfflength_mod.f
    fullmodel.so  # compiled f2py module
    fglists.dat

Python wrappers for code.  `models_all.py` allows fitting to both simple and
complex models, uses f2py to generate wrapper function for
`FullEfflength_mod.f` and replicates `Widthfun.py` fits.  `snr_catalog.py`
holds most physical parameters and synchrotron constants.

    models_all.py
    snr_catalog.py

Extras -- short FORTRAN code for learning, skeleton code for full port of
Sean's code over to Python (rather than using f2py wrapper), and a folder for
test runs of code (to be deleted).

    old-code/FortranTutorial.f
    old-code/FullEfflength_port.py
    test

Symlinks to matplotlibrc and fplot.py are just for plotting convenience.
A strange code bug is 


Files originally sent from Sean Ressler:
----------------------------------------
Main model code, model code with magnetic damping, simpler model code
(equation 6), and synchrotron emissivity table from Pacholczyk (197)

    src-doc-orig/FullEfflength.f
    src-doc-orig/FullefflengthPohlab.f
    src-doc-orig/Widthfun.py
    src-doc-orig/fglists.dat

Documentation of Sean's original code

    src-doc-orig/FWHMGuide.pdf
    src-doc-orig/sean-emails.md



