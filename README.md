README
======

Aaron Tran (supervisors: Rob Petre, Brian J. Williams)<br>
CRESST, NASA GSFC<br>

Analysis of thin synchrotron emitting filaments in historical supernova
remnants.  Pipeline to generate and fit radial intensity profiles from Chandra
observations of Tycho's SNR.  Model code (Sean Ressler and Steve Reynolds) to
compute magnetic fields and diffusion coefficients/scaling consistent with
observed rim width - energy relationships.

Work largely follows Ressler, Katsuda, Reynolds et al., ApJ, 2014

Repository
----------
Many necessary files (image files, raw Chandra data, energy spectra) are not
under version control and need to be redownloaded / regenerated to retrace or
extend this work.

Contact: forename.surname@berkeley.edu (substitute appropriately).

Most important folders are:

    code\
    code-profiles\
    code-models\
    data\

The file `pipeline.md` sketches out how to use the code pipeline to create
derived measurement products and run models for filament synchrotron emission.
