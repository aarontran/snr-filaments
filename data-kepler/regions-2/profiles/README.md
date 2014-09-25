README
======

Readme file for profiles generated 2014 September 25
Basically, the standard.  Using hanning window of length 13; binsize 1, etc.
as usual.  Region 1 looks to have too few counts!  Regions 8/9/10/11 are of
course sparse, but nothing really to be done there.

python ~/snr-research/code/profiles/prep_profile_fit.py ../regions-2.physreg prf prf-cts prf-proc -b 1 -l 0.7-1kev 1-2kev 2-4kev 4-7kev -w hanning -n 13 -v

Region 01: length=36.79px, thickness=14.00px
  Band 0.7-1kev: peak cts 91.0 in fit domain
  Band 1-2kev: peak cts 372.0 in fit domain
  Band 2-4kev: peak cts 103.0 in fit domain
  Band 4-7kev: peak cts 29.0 in fit domain
Region 02: length=37.11px, thickness=10.00px
  Band 0.7-1kev: peak cts 124.0 in fit domain
  Band 1-2kev: peak cts 506.0 in fit domain
  Band 2-4kev: peak cts 171.0 in fit domain
  Band 4-7kev: peak cts 58.0 in fit domain
Region 03: length=37.12px, thickness=8.00px
  Band 0.7-1kev: peak cts 93.0 in fit domain
  Band 1-2kev: peak cts 470.0 in fit domain
  Band 2-4kev: peak cts 207.0 in fit domain
  Band 4-7kev: peak cts 46.0 in fit domain
Region 04: length=33.92px, thickness=8.00px
  Band 0.7-1kev: peak cts 110.0 in fit domain
  Band 1-2kev: peak cts 495.0 in fit domain
  Band 2-4kev: peak cts 192.0 in fit domain
  Band 4-7kev: peak cts 44.0 in fit domain
Region 05: length=32.33px, thickness=8.00px
  Band 0.7-1kev: peak cts 101.0 in fit domain
  Band 1-2kev: peak cts 472.0 in fit domain
  Band 2-4kev: peak cts 172.0 in fit domain
  Band 4-7kev: peak cts 45.0 in fit domain
Region 06: length=28.76px, thickness=7.00px
  Band 0.7-1kev: peak cts 105.0 in fit domain
  Band 1-2kev: peak cts 420.0 in fit domain
  Band 2-4kev: peak cts 192.0 in fit domain
  Band 4-7kev: peak cts 54.0 in fit domain
Region 07: length=30.93px, thickness=7.00px
  Band 0.7-1kev: peak cts 77.0 in fit domain
  Band 1-2kev: peak cts 378.0 in fit domain
  Band 2-4kev: peak cts 158.0 in fit domain
  Band 4-7kev: peak cts 44.0 in fit domain
Region 08: length=29.39px, thickness=5.00px
  Band 0.7-1kev: peak cts 58.0 in fit domain
  Band 1-2kev: peak cts 245.0 in fit domain
  Band 2-4kev: peak cts 102.0 in fit domain
  Band 4-7kev: peak cts 29.0 in fit domain
Region 09: length=30.28px, thickness=12.00px
  Band 0.7-1kev: peak cts 85.0 in fit domain
  Band 1-2kev: peak cts 299.0 in fit domain
  Band 2-4kev: peak cts 122.0 in fit domain
  Band 4-7kev: peak cts 23.0 in fit domain
Region 10: length=30.61px, thickness=12.00px
  Band 0.7-1kev: peak cts 52.0 in fit domain
  Band 1-2kev: peak cts 242.0 in fit domain
  Band 2-4kev: peak cts 79.0 in fit domain
  Band 4-7kev: peak cts 20.0 in fit domain
Region 11: length=35.57px, thickness=11.00px
  Band 0.7-1kev: peak cts 77.0 in fit domain
  Band 1-2kev: peak cts 228.0 in fit domain
  Band 2-4kev: peak cts 88.0 in fit domain
  Band 4-7kev: peak cts 18.0 in fit domain
Wrote fit domain cuts to prf-proc_fit_cuts.dat, prf-proc_fit_cuts.npz
Done!
