Eyeballed fit parameters
========================
2014 December 9

Shape analysis "bounds" and eyeball fit parameters are obtained from
`proc_radio_xray_prfs.ipynb` in `code/models`, and should be fully documented
there as well.  Use that notebook to quickly replot the "fits" for any region
desired, and generate fit plots for the manuscript.

## (soft) "bounds" from shape analysis

Region 1 (A): minimum B0 = 20.910, maximum ab = 0.027
Region 2 (B): minimum B0 = 13.513, maximum ab = 0.042 (plateau)
Region 3 (C): minimum B0 = 13.513, maximum ab = 0.024 (plateau)
Region 4 (D): minimum B0 = 13.513, maximum ab = 0.042 (plateau)
Region 5 (E): minimum B0 = 20.290, maximum ab = 0.033
Region 6 (F): minimum B0 = 28.256, maximum ab = 0.019
Region 7 (G): minimum B0 = 19.688, maximum ab = 0.015
Region 8 (H): minimum B0 = 19.104, maximum ab = 0.027
Region 9 (I): minimum B0 = 24.676, maximum ab = 0.014
Region 10 (J): minimum B0 = 13.513, maximum ab = 0.059 (continuous rise)
Region 11 (K): minimum B0 = 31.396, maximum ab = 0.011
Region 12 (L): minimum ab ~ 0.1 (continuous rise)
Region 13 (M): minimum ab ~ 0.1 (continuous rise)
Region 14 (N): minimum B0 = 22.208, maximum ab = 0.012
Region 15 (O): minimum B0 = 20.598, maximum ab = 0.024
Region 16 (P): minimum B0 = 27.008, maximum ab = 0.010

## From manual (eyeballed) fitting procedure, with comments:

Region  1 (A): B0 =  50, ab = 0.02               DECENT rim example
Region  2 (B): B0 = 200, ab = 0.05     (plateau) GOOD plateau example
Region  3 (C): B0 =  15, ab = 0.01     (plateau) bad, few x-ray counts, hard to constrain, (ab, B0) could be (0.05, 60)
Region  4 (D): B0 = 100, ab = 0.05     (plateau) bad, same problem as above.  (ab, B0) could be (0.01, 15)
Region  5 (E): B0 = 120, ab = 0.03               bad, could be 2 filament structure
Region  6 (F): B0 = 300, ab = 0.025              bad, could be 2 filament structure
Region  7 (G): B0 = 250, ab = 0.02               DECENT, hard to fit radio and x-ray jointly
Region  8 (H): B0 = 300, ab = 0.02               DECENT, but x-ray rim hard to fit in general
Region  9 (I): B0 = 250, ab = 0.02               bad radio rim, fit doesn't work.
Region 10 (J): B0 = 300, ab = 0.3      (rise)    bad rise example -- best fit doesn't capture almost linear rise
Region 11 (K): B0 = 400, ab = 0.01               GOOD rim example
Region 12 (L): B0 = 250, ab = 0.2      (rise)    DECENT rise example
Region 13 (M): B0 = 200, ab = loss-lim (rise)    GOOD rise example
Region 14 (N): B0 =  23, ab = 0.004              bad, double filament on the NW
               B0 = 800, ab = 0.01                  alternate w/ extreme field gets both, but not plateau behind radio
Region 15 (O): B0 = 200, ab = 0.005              bad rim/plateau example (funny bump)
Region 16 (P): B0 = 150, ab = 0.012              bad, maybe multiple filaments (splaying apart)

I think Region B (plateau), K (rim), M (rise) are the best examples.
Region I is a nice example of where the radio rim does NOT work.
