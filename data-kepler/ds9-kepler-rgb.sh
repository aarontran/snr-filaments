#!/bin/bash
echo "Region file to load?"
read -e reg

ds9 -rgb \
    -red 0.7-1kev_mosaic.fits \
        -scale limits 4e-8 3.75e-6 \
    -green 1-2kev_mosaic.fits \
        -scale limits 2e-8 2.2e-6 \
    -blue 2-7kev_mosaic.fits \
        -scale limits 5e-9 4.5e-7 \
    -rgb lock scale yes \
    -log \
    -regions $reg
