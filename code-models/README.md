README
======
2014 July 26

Code is structured roughly as follows (from low-high level)

1. full model code, in Python and Fortran
2. `models.py` and `snr_catalog.py` wrap the fitting model / SNR details
3. `models_exec.py` calls the fits, builds SNRs
4. `models_disp.py` factory function for fit generation, and tables/plots
5. `*.ipynb` actually call functions with various data, settings, etc

To check/reproduce/rerun fits, simply go to the iPython notebooks.
For more involved things (running specific calculations, manipulating stuff to
display), move down the stack accordingly.

Tables for full model fits are in `tables/`, either scripts to regenerate
tables or the actual pkl files (if you got a copy from me directly, as the pkl
files are not version controlled).

To sanity check the fitting and error computation,
run `$ nosetests -v` from the terminal.  It's not comprehensive, but at least
verifies that (simple) fit results are consistent w/Sean's results, and that
our error procedure correctly gives +1 to chisqr


Main model code
---------------
Please use/read `FullEfflength_port.py`.  This is a Python port of Sean's model
code that should be more accurate; it does rely on the Fortran code
`FullEfflength_mod.f` for electron distribution calculation, which is the
slowest part of the code (run it through a profiler...).

For what the code is doing, refer to Sean's guide and my own code review(s)
(one level up).

    FullEfflength_port.py
    FullEfflength_mod.f
    fullmodel_recompile.py  # Code snippet to recompile FullEfflength_mod.f
    fullmodel.so  # compiled f2py module
    fglists.dat

Python wrappers for code.  `models.py` allows fitting to both simple and
complex models, uses f2py to generate wrapper function for
`FullEfflength_mod.f` and replicates `Widthfun.py` fits.  `snr_catalog.py`
holds most physical parameters and synchrotron constants.

    models.py
    snr_catalog.py

To get things done (perform fits, compute errors manually, etc):

    models_exec.py
    models_disp.py

Extras -- short FORTRAN code for learning, skeleton code for full port of
Sean's code over to Python (rather than using f2py wrapper), and a folder for
test runs of code (to be deleted).

    old-code/FortranTutorial.f
    old-code/FullEfflength_port.py
    test

Symlinks to matplotlibrc and fplot.py are just for plotting convenience.


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



