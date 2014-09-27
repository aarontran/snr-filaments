#!/bin/bash
echo "Region file to load?"
read -e reg

ds9 -rgb \
    -red 0.7-1kev_mosaic.fits \
        -scale limits 1.5e-8 1.9e-6 \
        -cmap value 1.95 0.22 \
    -green 1-2kev_mosaic.fits \
        -scale limits 8e-8 1e-5 \
        -cmap value 1.64 0.23 \
    -blue 2-7kev_mosaic.fits \
        -scale limits 4e-8 1e-5 \
        -cmap value 1.61 0.28 \
    -rgb lock scale yes \
    -log \
    -regions $reg
