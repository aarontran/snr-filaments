#!/bin/sh

punlearn specextract
#pset specextract grouptype=BIN_WIDTH binspec=.001
pset specextract weight=yes correct=no
pset specextract clobber=yes


specextract "/Users/ncoyle/Documents/ChandraData/data/10093/repro/acisf10093_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10093/repro/${1}10093.reg)]" ${1}10093

specextract "/Users/ncoyle/Documents/ChandraData/data/10094/repro/acisf10094_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10094/repro/${1}10094.reg)]" ${1}10094

specextract "/Users/ncoyle/Documents/ChandraData/data/10095/repro/acisf10095_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10095/repro/${1}10095.reg)]" ${1}10095

specextract "/Users/ncoyle/Documents/ChandraData/data/10096/repro/acisf10096_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10096/repro/${1}10096.reg)]" ${1}10096

specextract "/Users/ncoyle/Documents/ChandraData/data/10097/repro/acisf10097_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10097/repro/${1}10097.reg)]" ${1}10097

specextract "/Users/ncoyle/Documents/ChandraData/data/10902/repro/acisf10902_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10902/repro/${1}10902.reg)]" ${1}10902

specextract "/Users/ncoyle/Documents/ChandraData/data/10903/repro/acisf10903_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10903/repro/${1}10903.reg)]" ${1}10903

specextract "/Users/ncoyle/Documents/ChandraData/data/10904/repro/acisf10904_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10904/repro/${1}10904.reg)]" ${1}10904

specextract "/Users/ncoyle/Documents/ChandraData/data/10906/repro/acisf10906_repro_evt2.fits[sky=region(/Users/ncoyle/Documents/ChandraData/data/10906/repro/${1}10906.reg)]" ${1}10906


mathpha units=C expr="${1}10093.pi+${1}10094.pi+${1}10095.pi+${1}10096.pi+${1}10097.pi+${1}10902.pi+${1}10903.pi+${1}10904.pi+${1}10906.pi" outfil=${1}.pi exposure=CALC backscal='%'

addrmf "${1}10093.rmf ${1}10094.rmf ${1}10095.rmf ${1}10096.rmf ${1}10097.rmf ${1}10902.rmf ${1}10903.rmf ${1}10904.rmf ${1}10906.rmf" '0.1612156 0.1225565 0.2361635 0.1440111 0.1463404 0.053847516 0.03258367 0.04726812 0.056013404' ${1}.wrmf

addarf "${1}10093.arf ${1}10094.arf ${1}10095.arf ${1}10096.arf ${1}10097.arf ${1}10902.arf ${1}10903.arf ${1}10904.arf ${1}10906.arf" '0.1612156 0.1225565 0.2361635 0.1440111 0.1463404 0.053847516 0.03258367 0.04726812 0.056013404' ${1}.warf

#dmhedit infile=${1}.pi filelist="" operation=add key=BACKFILE value=annulus.pi
dmhedit infile=${1}.pi filelist="" operation=add key=RESPFILE value=${1}.wrmf
dmhedit infile=${1}.pi filelist="" operation=add key=ANCRFILE value=${1}.warf
