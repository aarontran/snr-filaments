Wednesday 2014 July 9
=====================
Short meeting with Brian in mid-late afternoon, to look over multipanel plots
of regions + spectra of subdivided regions.

Plots currently presented:
1. radial profiles + best fits for all bands
2. cuts for spectra extraction overlaid on 4.5-7keV band (highest energy)
3. spectra + residuals from downstream + rim regions (set by cut locations)

Discussion
----------
* Yeah, don't worry about the FWHM thing... (it overshoots, undershoots, can't
  really tell anyways, seriously)
* Spectrum normalization (using ct/s/keV/cm^2 instead of ct/s/keV) -- Brian
  hasn't seen that before, hasn't done the calculation himself.  It's usually
  just the XSPEC standard, "normalized counts /second /keV".  Ask Satoru about
  this later if important, can leave it be for now.
* Yes, get the figures/tables done first, before Sean's code.  We're at a point
  where we all need to be on the same page!


Things to adjust, presentation
------------------------------
* Change up colors slightly (done), for clarity.  Add shading to indicate where
  spectra are extracted.
* Remove "capped-FWHM" lines (done)
* we can put the plots on a website instead, doesn't have to be a PDF doc.
  Something for collaborators only, or what have you.  nbviewer looks okay.
* For spectra -- plot both rims + downstream areas with absorbed power law fits
  alone, no Gaussian lines.  Want to show that fits are worse for downstream
  region without fitting for thermal lines. (done)


Things to add / next steps
--------------------------
* Yes, go ahead and make your region changes / whatnot before sending.
* Add chi-squared values to spectra, on plot (preferable) or next to it.
* Make tickmarks readable, fiddle with formatting/sizing
* Calculate spectra with all 750 ks of data.  Ask Nina for her modified copy
  of Brian's script, which she got to work on a mac.
* Add tables of FWHM values, for people to see.  Eventually it will have to be
  in LaTeX.
