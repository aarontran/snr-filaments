#!/bin/sh

# Shell script to generate regions for Tycho's SNR
# Modified from script by Brian J. Williams
# Aaron Tran, 2014 July 12

# Call from data/regions-[n]/ directory

obsids=(10093 10094 10095 10096 10097 10902 10903 10904 10906)

mkdir tmp

for id in "${obsids[@]}";
do
    xpaset -p ds9 file "../../chandra/imgs_tycho/output${id}.fits"
    xpaset -p ds9 regions format ciao
    xpaset -p ds9 regions coord wcs
    xpaset -p ds9 regions load "${1}.reg"
    xpaset -p ds9 regions coord physical
    xpaset -p ds9 regions save "${1}_${id}.reg"
done
