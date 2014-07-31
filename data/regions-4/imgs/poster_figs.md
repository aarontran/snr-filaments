2014 July 29
============

Parameters used to generate poster figures (`tycho_*.tiff`)

## `4-7kev_mosaic.fits`

log scaling, cool colormap; scale limits 1.25e-8, 2.5e-7

    Try: 6.324, 64.143 (ra, dec in deg), 0.165 deg square --- same as the img generated for the
    SNR plot website (showing 3 color RGB image + regions).  This crop leaves space on top/right
    for the regions and labels.

    Use File/Export as .tiff, to get maximum resolution.

    2nd revision, brightening image (as suggested by Brian):
    Scale limits: 1.25e-8, 3e-07; contrast 1.3, bias 0.4; log scaling still
    Then crop and export as before.

    3rd revision... (2014 July 30)
    Contrast: 2.425
    Bias: 0.24548
    Smoothing: Gaussian smooth with radius 2
    Colormap: b
    Scaling: square root with limits 5e-9, 1.25e-7
    Crop: as before (centered on 6.324, 64.143; square 0.165 deg)

## 1.7 to 2.5 keV

    pset merge_obs infiles=10093/,10094/,10095/,10096/,10097/,10902/,10903/,10904/,10906/
    pset merge_obs outroot=test/thermal
    pset merge_obs bands=1.7:2.5:2.1
    pset merge_obs binsize=1

    Started at 00:54 (am).  Would be more efficient to split reproject/stack, but
    I don't anticipate running merge_obs again for Tycho data.
    Finished ~01:29.
    
    Colormap: bb, scale: square root, scale limits: 3e-08, 1.3e-06
    Crop: (ra, dec) = (6.324, 64.143) (degrees), 0.165 degrees square.
    File/Export to .tiff


## Large RGB image
My default setting is nice, but a little over exposed on the northwest.
Regions: use regions-4-img.reg (yellow color + width=2)... + region labels (times, fontsize 16, white)

    Here are the new parameters I'm using:
    
    0.7-1 keV limis: 2e-8 5e-7
    1-2 keV limits: 2e-8 5e-7
    2-7 keV limits: 2e-8 3e-7 (deliberately trying to bring out the thin rims)
    with asinh scaling

    Crop: (ra, dec) = (6.324, 64.143) (degrees), 0.165 degrees square.
    
    To account for poster scaling: multiply fontsizes by 2 (fontsize=32), set region widths = 3.
    
    NEW scalebar: width=4 line of approx 2.00 arcminutes.  2' fontsize=48
    
    Set frame parameters: 1208 x 1207 to match other images

## Large RGB image, round two
Borrowing settings from Nina's Tycho image:

    Channel Scale       Contrast    Bias
    Red     0, 1.5e-06  2.02        0.27
    Green   0, 2.2e-06  2.2         0.256
    Blue    0, 1.4e-06  2.7         0.19

    with linear scaling

    FINAL(?) settings -- had to twiddle scale params, units seem to be a bit
    different (dunno why)
    Channel     Scale           Contrast    Bias
    Red         1e-8, 4e-7      2.02        0.27
    Green       7e-9, 5e-7      2.2         0.256
    Blue        7e-9, 3.3e-7    2.7         0.19

    Apply same crop as before: 
    (ra, dec) = (6.324, 64.143) (degrees), 0.165 degrees square

    Export as tiff (no regions) and save image w/ 1208x1207 frame, zoom=1
    (and no colorbar, as before, mind you)
