README Cas A profiles
=====================

One-off calculation to get idea of counts/data.  We have oodles of counts, but
it doesn't really extend up to 9 keV.  But nice to see that we have excellent
counts to 4.5-7keV!

    python ~/snr-research/code/regions/ds9projplotter.py ../regions-0.reg prf -f
    ../../0.7-1kev_mosaic.fits ../../1-2kev_mosaic.fits ../../2-3kev_mosaic.fits
    ../../3-4.5kev_mosaic.fits ../../4.5-7kev_mosaic.fits ../../7-9kev_mosaic.fits
    -l 0.7-1kev 1-2kev 2-3kev 3-4.5kev 4.5-7kev 7-9kev -v -p plots/prf -s

    python ~/snr-research/code/regions/ds9projplotter.py ../regions-0.reg prf-cts
    -f ../../0.7-1kev_counts.fits ../../1-2kev_counts.fits ../../2-3kev_counts.fits
    ../../3-4.5kev_counts.fits ../../4.5-7kev_counts.fits ../../7-9kev_counts.fits
    -l 0.7-1kev 1-2kev 2-3kev 3-4.5kev 4.5-7kev 7-9kev -v

    python ~/snr-research/code/profiles/prep_profile_fit.py \
           ../regions-0.physreg prf prf-cts prf-proc --binsize 1 \
           --labels 0.7-1kev 1-2kev 2-3kev 3-4.5kev 4.5-7kev 7-9kev \
           --window hanning --window-n 11 --verbose


Region 01: length=30.62px, thickness=7.00px
  Band 0.7-1kev: peak cts 112.0 in fit domain
  Band 1-2kev: peak cts 904.0 in fit domain
  Band 2-3kev: peak cts 269.0 in fit domain
  Band 3-4.5kev: peak cts 156.0 in fit domain
  Band 4.5-7kev: peak cts 50.0 in fit domain
  Band 7-9kev: peak cts 6.0 in fit domain
Region 02: length=31.56px, thickness=6.00px
  Band 0.7-1kev: peak cts 87.0 in fit domain
  Band 1-2kev: peak cts 895.0 in fit domain
  Band 2-3kev: peak cts 301.0 in fit domain
  Band 3-4.5kev: peak cts 175.0 in fit domain
  Band 4.5-7kev: peak cts 63.0 in fit domain
  Band 7-9kev: peak cts 11.0 in fit domain
Region 03: length=31.19px, thickness=6.00px
  Band 0.7-1kev: peak cts 73.0 in fit domain
  Band 1-2kev: peak cts 723.0 in fit domain
  Band 2-3kev: peak cts 263.0 in fit domain
  Band 3-4.5kev: peak cts 157.0 in fit domain
  Band 4.5-7kev: peak cts 73.0 in fit domain
  Band 7-9kev: peak cts 7.0 in fit domain
Region 04: length=31.04px, thickness=6.00px
  Band 0.7-1kev: peak cts 76.0 in fit domain
  Band 1-2kev: peak cts 693.0 in fit domain
  Band 2-3kev: peak cts 278.0 in fit domain
  Band 3-4.5kev: peak cts 163.0 in fit domain
  Band 4.5-7kev: peak cts 66.0 in fit domain
  Band 7-9kev: peak cts 8.0 in fit domain
Region 05: length=31.94px, thickness=6.00px
  Band 0.7-1kev: peak cts 75.0 in fit domain
  Band 1-2kev: peak cts 712.0 in fit domain
  Band 2-3kev: peak cts 236.0 in fit domain
  Band 3-4.5kev: peak cts 166.0 in fit domain
  Band 4.5-7kev: peak cts 75.0 in fit domain
  Band 7-9kev: peak cts 6.0 in fit domain
Region 06: length=31.04px, thickness=6.00px
  Band 0.7-1kev: peak cts 79.0 in fit domain
  Band 1-2kev: peak cts 753.0 in fit domain
  Band 2-3kev: peak cts 281.0 in fit domain
  Band 3-4.5kev: peak cts 161.0 in fit domain
  Band 4.5-7kev: peak cts 62.0 in fit domain
  Band 7-9kev: peak cts 7.0 in fit domain
Region 07: length=52.41px, thickness=7.00px
  Band 0.7-1kev: peak cts 71.0 in fit domain
  Band 1-2kev: peak cts 715.0 in fit domain
  Band 2-3kev: peak cts 316.0 in fit domain
  Band 3-4.5kev: peak cts 206.0 in fit domain
  Band 4.5-7kev: peak cts 104.0 in fit domain
  Band 7-9kev: peak cts 13.0 in fit domain
Region 08: length=42.88px, thickness=7.00px
  Band 0.7-1kev: peak cts 65.0 in fit domain
  Band 1-2kev: peak cts 824.0 in fit domain
  Band 2-3kev: peak cts 383.0 in fit domain
  Band 3-4.5kev: peak cts 330.0 in fit domain
  Band 4.5-7kev: peak cts 143.0 in fit domain
  Band 7-9kev: peak cts 17.0 in fit domain
Region 09: length=42.26px, thickness=7.00px
  Band 0.7-1kev: peak cts 58.0 in fit domain
  Band 1-2kev: peak cts 842.0 in fit domain
  Band 2-3kev: peak cts 379.0 in fit domain
  Band 3-4.5kev: peak cts 329.0 in fit domain
  Band 4.5-7kev: peak cts 133.0 in fit domain
  Band 7-9kev: peak cts 15.0 in fit domain
Region 10: length=40.15px, thickness=8.00px
  Band 0.7-1kev: peak cts 69.0 in fit domain
  Band 1-2kev: peak cts 765.0 in fit domain
  Band 2-3kev: peak cts 358.0 in fit domain
  Band 3-4.5kev: peak cts 287.0 in fit domain
  Band 4.5-7kev: peak cts 137.0 in fit domain
  Band 7-9kev: peak cts 14.0 in fit domain
Region 11: length=78.92px, thickness=21.00px
  Band 0.7-1kev: peak cts 92.0 in fit domain
  Band 1-2kev: peak cts 853.0 in fit domain
  Band 2-3kev: peak cts 265.0 in fit domain
  Band 3-4.5kev: peak cts 163.0 in fit domain
  Band 4.5-7kev: peak cts 64.0 in fit domain
  Band 7-9kev: peak cts 10.0 in fit domain
Region 12: length=54.08px, thickness=12.00px
  Band 0.7-1kev: peak cts 55.0 in fit domain
  Band 1-2kev: peak cts 574.0 in fit domain
  Band 2-3kev: peak cts 328.0 in fit domain
  Band 3-4.5kev: peak cts 259.0 in fit domain
  Band 4.5-7kev: peak cts 127.0 in fit domain
  Band 7-9kev: peak cts 12.0 in fit domain
Region 13: length=65.68px, thickness=30.00px
  Band 0.7-1kev: peak cts 142.0 in fit domain
  Band 1-2kev: peak cts 1911.0 in fit domain
  Band 2-3kev: peak cts 684.0 in fit domain
  Band 3-4.5kev: peak cts 249.0 in fit domain
  Band 4.5-7kev: peak cts 73.0 in fit domain
  Band 7-9kev: peak cts 15.0 in fit domain
Region 14: length=68.75px, thickness=15.00px
  Band 0.7-1kev: peak cts 178.0 in fit domain
  Band 1-2kev: peak cts 4178.0 in fit domain
  Band 2-3kev: peak cts 2220.0 in fit domain
  Band 3-4.5kev: peak cts 1679.0 in fit domain
  Band 4.5-7kev: peak cts 646.0 in fit domain
  Band 7-9kev: peak cts 34.0 in fit domain
Region 15: length=46.17px, thickness=12.00px
  Band 0.7-1kev: peak cts 255.0 in fit domain
  Band 1-2kev: peak cts 2379.0 in fit domain
  Band 2-3kev: peak cts 826.0 in fit domain
  Band 3-4.5kev: peak cts 517.0 in fit domain
  Band 4.5-7kev: peak cts 208.0 in fit domain
  Band 7-9kev: peak cts 19.0 in fit domain
Region 16: length=36.57px, thickness=12.00px
  Band 0.7-1kev: peak cts 152.0 in fit domain
  Band 1-2kev: peak cts 1667.0 in fit domain
  Band 2-3kev: peak cts 657.0 in fit domain
  Band 3-4.5kev: peak cts 577.0 in fit domain
  Band 4.5-7kev: peak cts 255.0 in fit domain
  Band 7-9kev: peak cts 19.0 in fit domain
Region 17: length=29.07px, thickness=12.00px
  Band 0.7-1kev: peak cts 334.0 in fit domain
  Band 1-2kev: peak cts 3944.0 in fit domain
  Band 2-3kev: peak cts 1467.0 in fit domain
  Band 3-4.5kev: peak cts 1060.0 in fit domain
  Band 4.5-7kev: peak cts 395.0 in fit domain
  Band 7-9kev: peak cts 23.0 in fit domain
Region 18: length=40.03px, thickness=12.00px
  Band 0.7-1kev: peak cts 158.0 in fit domain
  Band 1-2kev: peak cts 1018.0 in fit domain
  Band 2-3kev: peak cts 268.0 in fit domain
  Band 3-4.5kev: peak cts 156.0 in fit domain
  Band 4.5-7kev: peak cts 53.0 in fit domain
  Band 7-9kev: peak cts 7.0 in fit domain
Region 19: length=28.46px, thickness=10.00px
  Band 0.7-1kev: peak cts 232.0 in fit domain
  Band 1-2kev: peak cts 1685.0 in fit domain
  Band 2-3kev: peak cts 631.0 in fit domain
  Band 3-4.5kev: peak cts 426.0 in fit domain
  Band 4.5-7kev: peak cts 188.0 in fit domain
  Band 7-9kev: peak cts 15.0 in fit domain
Region 20: length=24.26px, thickness=15.00px
  Band 0.7-1kev: peak cts 723.0 in fit domain
  Band 1-2kev: peak cts 6600.0 in fit domain
  Band 2-3kev: peak cts 2181.0 in fit domain
  Band 3-4.5kev: peak cts 1314.0 in fit domain
  Band 4.5-7kev: peak cts 449.0 in fit domain
  Band 7-9kev: peak cts 23.0 in fit domain
Region 21: length=37.74px, thickness=10.00px
  Band 0.7-1kev: peak cts 318.0 in fit domain
  Band 1-2kev: peak cts 3358.0 in fit domain
  Band 2-3kev: peak cts 1299.0 in fit domain
  Band 3-4.5kev: peak cts 824.0 in fit domain
  Band 4.5-7kev: peak cts 321.0 in fit domain
  Band 7-9kev: peak cts 26.0 in fit domain
Region 22: length=35.45px, thickness=9.00px
  Band 0.7-1kev: peak cts 186.0 in fit domain
  Band 1-2kev: peak cts 2318.0 in fit domain
  Band 2-3kev: peak cts 1059.0 in fit domain
  Band 3-4.5kev: peak cts 900.0 in fit domain
  Band 4.5-7kev: peak cts 442.0 in fit domain
  Band 7-9kev: peak cts 39.0 in fit domain
Region 23: length=21.68px, thickness=10.00px
  Band 0.7-1kev: peak cts 331.0 in fit domain
  Band 1-2kev: peak cts 3333.0 in fit domain
  Band 2-3kev: peak cts 1329.0 in fit domain
  Band 3-4.5kev: peak cts 874.0 in fit domain
  Band 4.5-7kev: peak cts 312.0 in fit domain
  Band 7-9kev: peak cts 25.0 in fit domain
Region 24: length=44.33px, thickness=7.00px
  Band 0.7-1kev: peak cts 245.0 in fit domain
  Band 1-2kev: peak cts 3228.0 in fit domain
  Band 2-3kev: peak cts 1472.0 in fit domain
  Band 3-4.5kev: peak cts 1077.0 in fit domain
  Band 4.5-7kev: peak cts 490.0 in fit domain
  Band 7-9kev: peak cts 41.0 in fit domain
Region 25: length=35.14px, thickness=6.00px
  Band 0.7-1kev: peak cts 171.0 in fit domain
  Band 1-2kev: peak cts 1930.0 in fit domain
  Band 2-3kev: peak cts 910.0 in fit domain
  Band 3-4.5kev: peak cts 699.0 in fit domain
  Band 4.5-7kev: peak cts 302.0 in fit domain
  Band 7-9kev: peak cts 23.0 in fit domain
Region 26: length=32.20px, thickness=10.00px
  Band 0.7-1kev: peak cts 243.0 in fit domain
  Band 1-2kev: peak cts 1780.0 in fit domain
  Band 2-3kev: peak cts 619.0 in fit domain
  Band 3-4.5kev: peak cts 369.0 in fit domain
  Band 4.5-7kev: peak cts 142.0 in fit domain
  Band 7-9kev: peak cts 12.0 in fit domain
Region 27: length=23.70px, thickness=20.00px
  Band 0.7-1kev: peak cts 368.0 in fit domain
  Band 1-2kev: peak cts 3211.0 in fit domain
  Band 2-3kev: peak cts 1047.0 in fit domain
  Band 3-4.5kev: peak cts 709.0 in fit domain
  Band 4.5-7kev: peak cts 302.0 in fit domain
  Band 7-9kev: peak cts 27.0 in fit domain
Wrote fit domain cuts to prf-proc_fit_cuts.dat, prf-proc_fit_cuts.npz
Done!

