#!/bin/sh

# Shell script to generate regions for Tycho's SNR
# Modified from script by Brian J. Williams
# Aaron Tran, 2014 July 12

# First call heainit, then ds9&, then
# call script from data/regions-[n]/ directory

# obsids=(10093 10094 10095 10096 10097 10902 10903 10904 10906)
obsids=(6714 6715 6716 6717 6718 7366)

for id in "${obsids[@]}";
do
#    xpaset -p ds9 file "../../chandra/imgs_tycho/output${id}.fits"
    xpaset -p ds9 file "../../chandra/${id}/repro/acisf0${id}_repro_evt2.fits"
#    xpaset -p ds9 regions format ciao
    xpaset -p ds9 regions format ds9
    xpaset -p ds9 regions system wcs
    xpaset -p ds9 regions load "${1}.reg"
    xpaset -p ds9 regions format ciao
    xpaset -p ds9 regions system physical
    xpaset -p ds9 regions save "${1}_${id}.reg"
done
