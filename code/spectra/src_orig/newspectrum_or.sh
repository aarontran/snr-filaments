#!/bin/bash
#
#  Shell script to generate spectra for Tycho's SNR
#

punlearn specextract
pset specextract infile=@in${1}.lis
pset specextract outroot=@out${1}.lis
pset specextract grouptype=NONE
#pset specextract pbkfile=@pbkfile.lis
pset specextract asp=@aspfile.lis
pset specextract msk=@mskfile.lis
pset specextract badpixfile=@bpixfile.lis
pset specextract weight=yes correct=no
pset specextract grouptype=NONE binspec=NONE
pset specextract clobber=yes

specextract verbose=2

mathpha units=C expr="${1}10093.pi+${1}10094.pi+${1}10095.pi+${1}10096.pi+${1}10097.pi+${1}10902.pi+${1}10903.pi+${1}10904.pi+${1}10906.pi" outfil=${1}.pi exposure=CALC backscal='%'

addrmf "${1}10093.rmf ${1}10094.rmf ${1}10095.rmf ${1}10096.rmf ${1}10097.rmf ${1}10902.rmf ${1}10903.rmf ${1}10904.rmf ${1}10906.rmf"  '0.1612156 0.1225565 0.2361635 0.1440111 0.1463404 0.053847516 0.03258367 0.04726812 0.056013404' ${1}.rmf 

addarf "${1}10093.arf ${1}10094.arf ${1}10095.arf ${1}10096.arf ${1}10097.arf ${1}10902.arf ${1}10903.arf ${1}10904.arf ${1}10906.arf"  '0.1612156 0.1225565 0.2361635 0.1440111 0.1463404 0.053847516 0.03258367 0.04726812 0.056013404' ${1}.arf 

#dmhedit infile=${1}.pi filelist="" operation=add key=BACKFILE value=tychobg.pi
dmhedit infile=${1}.pi filelist="" operation=add key=RESPFILE value=${1}.wrmf 
dmhedit infile=${1}.pi filelist="" operation=add key=ANCRFILE value=${1}.warf 

