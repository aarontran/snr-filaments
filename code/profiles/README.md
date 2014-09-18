README for profile fitting code
===============================

Fitting functions, utilities
----------------------------
Add/modify the ipython notebook, then export to `.py` file
Reason being that the ipynb allows rich markdown w/ LaTeX, easier to understand
the fit functions being used.

    fit_funcs.ipynb
    fit_funcs.py

Nosetests in IPython

    ipython_nose-master
    ipython_nose.py

Fitting notebooks etc..
-----------------------

    fwhm-process.ipynb
    fwhm-sn1006.ipynb
    plotter_plots
    plotter_prfs_specs.ipynb
    plotter_prfs_specs_2014-07-14_email.html
    plotter_prfs_specs_2014-07-14_email.ipynb
    profile_process.ipynb
    profile_process_SN1006_2014_06_30.ipynb

Inside notebooks specify filenames.
They will spew out plots etc.  You twiddle the settings until you get nice
plots and FWHMs.

Then go ahead and save them to pickle files.

Utilities etc for fitting notebooks
-----------------------------------

    fsmooth.py
    crvfit.py
    ffwhm.py
