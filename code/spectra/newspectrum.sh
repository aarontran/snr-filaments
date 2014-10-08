#!/bin/sh

# Script to generate spectra for individual regions, for each ObsID
# Modified from script by Brian J. Williams, Nina Coyle
# Aaron Tran, 2014 July 12

# Run from same directory as regions from newregion.sh
# Input: filestem of original region file -- same as filestem of output spectra

punlearn specextract
pset specextract weight=yes correct=no
pset specextract grouptype=NONE binspec=NONE
pset specextract clobber=yes

# obsids=(10093 10094 10095 10096 10097 10902 10903 10904 10906)
# REMOVE/ADD LEADING ZERO TO acisf0${id}_repro_evt2.fits ...
# if using Tycho vs. Kepler/CasA data
obsids=(6714 6715 6716 6717 6718 7366)

# Specextract - this takes 2-4 hours
for id in "${obsids[@]}";
do
    specextract "../../chandra/${id}/repro/acisf0${id}_repro_evt2.fits[sky=@${1}_${id}.reg]" "spectra/${1}_${id}"
done
