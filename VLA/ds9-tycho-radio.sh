#!/bin/bash
echo "Region file to load?"
read -e reg

ds9 -view layout vertical \
    -colorbar orientation vertical \
    -rgb \
    -red TYCHO_IF1.FITS \
        -scale limits 0.00015 0.0035 \
        -cmap value 2.25 0.55 \
    -green ../data-tycho/1-1.7kev_mosaic.fits \
        -scale limits 7e-9 4e-7 \
        -cmap value 3.59 0.12 \
    -blue ../data-tycho/4-7kev_mosaic.fits \
        -scale limits 7e-9 4e-7 \
        -cmap value 3.59 0.12 \
    -rgb lock scale yes \
    -asinh \
    -regions $reg
