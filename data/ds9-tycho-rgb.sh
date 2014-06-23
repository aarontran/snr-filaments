#!/bin/bash
echo "Region file to load?"
read reg

ds9 -view layout vertical \
    -colorbar orientation vertical \
    -rgb \
    -red 0.7-1kev_mosaic.fits \
        -scale limits 1e-8 2.1e-6 \
    -green 1-2kev_mosaic.fits \
        -scale limits 7e-9 7e-7 \
    -blue 2-7kev_mosaic.fits \
        -scale limits 7e-9 5e-7 \
    -rgb lock scale yes \
    -histequ \
    -regions $reg
