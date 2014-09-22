#!/bin/bash
echo "Region file to load?"
read -e reg

ds9 -rgb \
    -red 0.7-1kev_mosaic.fits \
        -scale limits 2.5e-8 1.9e-6 \
    -green 1-2kev_mosaic.fits \
        -scale limits 1e-7 1e-5 \
    -blue 2-7kev_mosaic.fits \
        -scale limits 1e-7 1e-5 \
    -rgb lock scale yes \
    -log \
    -regions $reg
