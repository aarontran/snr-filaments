#!/bin/bash
echo "Region file to load?"
read reg

ds9 -view layout vertical \
    -colorbar orientation vertical \
    -rgb \
    -red 0.7-1kev_mosaic_bin.fits \
        -scale limits 4e-8 2.1e-6 \
    -green 1-2kev_mosaic_bin.fits \
        -scale limits 4e-8 2.42e-6 \
    -blue 2-7kev_mosaic_bin.fits \
        -scale limits 3e-8 1.5e-6 \
    -rgb lock scale yes \
    -histequ \
    -regions $reg
