"""
Short script to generate Tycho image in APLpy
Should not be run stand-alone (won't do much); run with:
$ ipython --pylab
and make any desired customizations.  Then, go ahead and save...

Also misc utilities/tidbits in functions here

Aaron Tran
"""

import aplpy

def main():
    """Load Tycho image and tweak settings for one-column manuscript fig"""

    f = aplpy.FITSFigure('ds9rgb.fits', figsize=(5,4.5))  # Saved from DS9
    f.show_rgb('ds9rgb.png')            # Saved from DS9
    f.show_regions('../regions-6/regions-6-box.reg')  # f.remove_layer(...)

    f.recenter(ra2deg(0,25,16.5), 64.143, width=10./60, height=10./60)
    # Slightly off-center so labels on N/W have more space (literally...)

    f.tick_labels.set_yformat('dd:mm')
    f.tick_labels.set_xformat('hh:mm:ss')
    f.axis_labels.set_ypad(-12)
    f.refresh()

    f_mpl = f.image.figure  # Seems to go straight to pyplot api
    f_mpl.tight_layout()

    return f


def clear_layers(f):
    """Clear all layers from an aplpy.FITSFigure"""
    k = f._layers.keys()  # Any better way to do this?
    map(f.remove_layer, k)

def ra2deg(h,m,s):
    """Convert RA specified in hh:mm:ss to degrees"""
    return (h+m/60.+s/3600.)/24.*360.


if __name__ == '__main__':
    main()
