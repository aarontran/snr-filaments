#!/bin/bash
echo "Region file to load?"
read -e reg

ds9 -view layout vertical \
    -colorbar orientation vertical \
    -fits TYCHO_IF1.FITS \
    -scale limits 0 0.0037 \
    -cmap bb \
    -cmap value 1.804 0.263 \
    -linear \
    -regions $reg
