#!/bin/bash

# Short script to batch convert ps2pdf and remove intermediate .ps files
# Aaron Tran, 2014 June 17

echo "Enter file stem (e.g., stem for ./plots/plt_src1_grp.ps is /plots/plt):"
read stem

nfiles=`ls -l ${stem}_src*_grp.ps | wc -l`
echo "${nfiles} .ps files identified"

for i in `seq 1 ${nfiles}`;
do
    infile=${stem}_src${i}_grp.ps
    outfile=${stem}_src${i}_grp.pdf
    echo "Converting ${infile} to ${outfile}"
    ps2pdf $infile $outfile
    rm $infile
done
echo Done!
