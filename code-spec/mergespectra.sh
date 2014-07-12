#!/bin/bash

# Script to merge output spectra from newregion.sh, newspectrum.sh
# Modified from script by Brian J. Williams, Nina Coyle
# Aaron Tran, 2014 July 12

# Input: filestem of spectra -- call script from same directory as spectra
# e.g., stem for regions4down10093.pi is 'regions4down'

# Check cmd line args
if [ "$#" -ne 1 ]
then
    echo "Error: one argument required, exiting"
    exit 1
fi

# Tycho's SNR, Hughes 2007
obsids=(10093 10094 10095 10096 10097 10902 10903 10904 10906)
wghts=(0.1612156 0.1225565 0.2361635 0.1440111 0.1463404 \
       0.053847516 0.03258367 0.04726812 0.056013404)

# Check number of regions to process
nregs=`ls -l "${1}_${obsids[0]}_src*.pi" | wc -l`
if [ ${nregs} != 0 ]
then
    echo "${nregs} regions' spectra will be merged"
else
    echo "No spectra to merge, exiting"
    exit 0
fi

# Perform merge operation for each region
for n in `seq 1 ${nregs}`;
do
    # Build argument strings -- list of source files
    pi_str=""
    rmf_str=""
    arf_str=""
    for id in "${obsids[@]}";
    do
        pi_str="${pi_str}+${1}_${id}_src${n}.pi"
        rmf_str="${rmf_str} ${1}_${id}_src${n}.rmf"
        arf_str="${arf_str} ${1}_${id}_src${n}.arf"
    done
    # Remove extra '+' operator
    pi_str=${pi_str:1}

    # Merge pi, rmf, arf files (FTOOLS)
    mathpha units=C expr=${pi_str} outfil=${1}_src${n}.pi exposure=CALC \
            backscal='%' areascal='%' ncomments=0
    addrmf "${rmf_str}" "${wghts[*]}" "${1}_src${n}.wrmf"
    addarf "${arf_str}" "${wghts[*]}" "${1}_src${n}.warf"

    # Link files (CIAO)
    dmhedit infile="${1}_src${n}.pi" filelist="" operation=add key=RESPFILE \
            value="${1}_src${n}.wrmf"
    dmhedit infile="${1}_src${n}.pi" filelist="" operation=add key=ANCRFILE \
            value="${1}_src${n}.warf"

    # Generate grouped spectra
    grppha infile="${1}_src${n}.pi" outfile="${1}_src${n}_grp.pi" comm="group min 15 & exit"
done

# Move the merged spectra up as you desire???..
# Here it may be good to clean up the spectra, then cd upwards and rm regions
# But give user choice/flag... incase merging needs to be vetted/redone
