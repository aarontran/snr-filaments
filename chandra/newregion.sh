#!/bin/sh
#
#  Shell script to generate regions for Tycho's SNR
#
xpaset -p ds9 file 10093/repro/acisf10093_repro_evt2.fits  # output10093.fits
xpaset -p ds9 regions format ciao           # acisf10096_repro_evt2.fits
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10093/repro/${1}10093.reg
#
xpaset -p ds9 file 10094/repro/acisf10094_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10094/repro/${1}10094.reg
#
xpaset -p ds9 file 10095/repro/acisf10095_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10095/repro/${1}10095.reg
#
xpaset -p ds9 file 10096/repro/acisf10096_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10096/repro/${1}10096.reg
#
xpaset -p ds9 file 10097/repro/acisf10097_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10097/repro/${1}10097.reg
#
xpaset -p ds9 file 10902/repro/acisf10902_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10902/repro/${1}10902.reg
#
xpaset -p ds9 file 10903/repro/acisf10903_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10903/repro/${1}10903.reg
#
xpaset -p ds9 file 10904/repro/acisf10904_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10904/repro/${1}10904.reg
#
xpaset -p ds9 file 10906/repro/acisf10906_repro_evt2.fits
xpaset -p ds9 regions format ciao
xpaset -p ds9 regions coord wcs
xpaset -p ds9 regions load ${1}coord.reg
xpaset -p ds9 regions coord physical
xpaset -p ds9 regions save 10906/repro/${1}10906.reg

sed s/sample/${1}/ insample.lis > spectra/in${1}.lis
#this just generates a list of prefixes for the output files
sed s/sample/${1}/ outsample.lis > spectra/out${1}.lis
