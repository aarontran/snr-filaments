"""
Short script to generate Tycho images in APLpy
Should not be run stand-alone (won't do much); run with:
$ ipython --pylab
and make any desired customizations.  Then, go ahead and save...

Aaron Tran
"""

import aplpy

def main():
    #manuscript_fig_snr()
    poster_fig_snr()

# ================================
# APLpy tidbits for poster figures
# ================================

#def poster_fig_snr_radio():
#    """Tycho 1.375 GHz image with scalebar, attempt to match size/scale of
#    other poster images"""
#    f = aplpy.FITSFigure('../TYCHO_IF1.FITS', figsize=(10, 10))
#    f.show_colorscale(vmin=0, vmax=0.0037, stretch='linear', cmap='hot')
#    apply_tycho_frame(f, epoch='B1950')
#
#    prep_save_pdf(f, 'snr_radio_poster.pdf')
#
#
#def poster_fig_snr_nonthermal():
#    """Tycho 4-7keV mosaic with scalebar, showing nonthermal rims"""
#    f = aplpy.FITSFigure('../4-7kev_mosaic.fits', figsize=(10, 10))
#    f.show_colorscale(vmin=5e-9, vmax=1.25e-7, stretch='sqrt', cmap='hot')
#    apply_tycho_frame(f, epoch='J2000')
#
#    prep_save_pdf(f, 'snr_4-7kev_poster.pdf')


def poster_fig_snr():
    """RGB figure of Tycho with X-ray region overlay, scalebar
    Uses regions w/ larger text, thicker boxes"""
    f = aplpy.FITSFigure('ds9rgb.fits', figsize=(10, 10))
    f.show_rgb('ds9rgb.png')
    f.show_regions('regs/regions-6-box-poster.reg')
    apply_tycho_frame(f, epoch='J2000')

    prep_save_pdf(f, 'snr_poster.pdf')


def prep_save_pdf(f, fname):
    """Image prepwork for poster figures, apply to all for consistency"""
    # Remove all ticks and associated labels
    f.ticks.hide()
    f.hide_tick_labels()
    f.hide_xaxis_label()
    f.hide_yaxis_label()

    f.image.figure.tight_layout()  # Seems to go straight to pyplot api

    # Add white 2 arcmin. scalebar
    f.add_scalebar(2./60)
    f.scalebar.set_label(r'$2^\prime$')
    f.scalebar.set_color('white')
    f.scalebar.set_linewidth(2)

    f.refresh()
    f.save(fname)  # Saves at native resolution of png/fits?


def apply_tycho_frame(f, epoch='J2000'):
    """Set a default centering and frame size to capture SNR + regions
    Slightly off-center so labels on N/W have more space (literally...)
    """
    if epoch == 'J2000':
        f.recenter(ra2deg(0,25,17.5), 64.141, width=10.5/60, height=10.5/60)
    elif epoch == 'B1950':
        # Convert J2000 to B1950 using
        # http://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm
        # Input: RA = 6.322917, DEC = 64.141
        # Output: RA = 5.623, DEC = 63.8641 (B1950)
        f.recenter(5.623, 63.8641, width=10.5/60, height=10.5/60)


# ====================================
# APLpy tidbits for manuscript figures
# ====================================

def manuscript_fig_snr():
    """One-column RGB figure of Tycho with X-ray regions-6 overlay"""
    f = aplpy.FITSFigure('ds9rgb.fits', figsize=(5,4.5))
    f.show_rgb('ds9rgb.png')
    f.show_regions('regs/regions-6-box.reg')

    f.recenter(ra2deg(0,25,16.5), 64.143, width=10./60, height=10./60)
    # Slightly off-center so labels on N/W have more space (literally...)

    f.tick_labels.set_yformat('dd:mm')
    f.tick_labels.set_xformat('hh:mm:ss')
    f.axis_labels.set_ypad(-12)
    f.refresh()

    f_mpl = f.image.figure  # Seems to go straight to pyplot api
    f_mpl.tight_layout()

    f.save('snr_manuscript.pdf', dpi=300)

    return f


# ===============
# Misc. utilities
# ===============

def clear_layers(f):
    """Clear all layers from an aplpy.FITSFigure"""
    k = f._layers.keys()  # Any better way to do this?
    map(f.remove_layer, k)

def ra2deg(h,m,s):
    """Convert RA specified in hh:mm:ss to degrees"""
    return (h+m/60.+s/3600.)/24.*360.


if __name__ == '__main__':
    main()
