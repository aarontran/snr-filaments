Code output post-mortem
=======================

Let's walk through what happened.
I ran it for ~8 hours last night and it only got through 2 filaments...
(error post mortem addressed on Friday Aug. 29)

Filament 1 tables
-----------------

mu	    eta2	                            B0	                            chisqr  NOTES
0.33	5.678 +77.19/-4.64 (std: ± 13.167)	127 +358/-27 (std: ± 52.660)	0.0568  eta2 upper lim. too low
0.50	3.849 +53.05/-2.88 (std: ± 6.242)	116 +658/-17.9 (std: ± 29.368)	0.0741  B0/eta2 upper lims. too low
1.00	2.590 +5.39/-1.68 (std: ± 2.598)	102 +17.6/-8.84 (std: ± 11.628)	0.1415
1.50	2.457 +3.21/-1.50 (std: ± 2.065)	97.1 +9.18/-5.87 (std: ± 7.133)	0.2221
2.00	2.640 +2.88/-1.57 (std: ± 2.067)	94.2 +6.26/-4.51 (std: ± 5.231)	0.3021

    \begin{tabular}{@{}rllr@{}}
    \toprule
    {} & \multicolumn{3}{c}{Filament 1} \\
    \cmidrule(l){2-4}
    $\mu$ (-) & $\eta_2$ (-) & $B_0$ ($\mu$G) & $\chi^2$\\
    \midrule
    0.33 & ${5.7}^{+77}_{-4.6}$ & ${127}^{+3.6 \times 10^{2}}_{-27}$ & 0.0568\\
    0.50 & ${3.8}^{+53}_{-2.9}$ & ${116}^{+6.6 \times 10^{2}}_{-18}$ & 0.0741\\
    1.00 & ${2.6}^{+5.4}_{-1.7}$ & ${102.5}^{+18}_{-8.8}$ & 0.1415\\
    1.50 & ${2.5}^{+3.2}_{-1.5}$ & ${97.1}^{+9.2}_{-5.9}$ & 0.2221\\
    2.00 & ${2.6}^{+2.9}_{-1.6}$ & ${94.2}^{+6.3}_{-4.5}$ & 0.3021\\
    \bottomrule

Filament 2 tables
-----------------

mu	    eta2	                            B0	                                chisqr      notes
0.33	0.015 +0.21/--0.21 (std: ± 0.330)	132 +2.33/-3.75 (std: ± 3.412)	    52.8783     B0 lower limit bad, prob. too small.  eta2 lower limit is inane
0.50	0.012 +0.36/-0.01 (std: ± 1.626)	131 +4.33/-0.393 (std: ± 35.371)	53.6116     Best fit values suspect? eta2 hit 0. for lower limit
1.00	0.176 +0.25/-0.18 (std: ± 0.263)	131 +2.19/-1.02 (std: ± 1.898)	    52.3628     eta2 hit 0. for lower limit
1.50	0.316 +0.16/-0.12 (std: ± 0.328)	130 +1.01/-0.747 (std: ± 1.897)	    51.0397     B0 lower limit suspect
2.00	0.343 +0.61/-0.31 (std: ± 0.361)	128 +2.55/-1.39 (std: ± 1.457)	    51.5343     B0 upper limit, slightly too small

    \begin{tabular}{@{}rllr@{}}
    \toprule
    {} & \multicolumn{3}{c}{Filament 2} \\
    \cmidrule(l){2-4}
    $\mu$ (-) & $\eta_2$ (-) & $B_0$ ($\mu$G) & $\chi^2$\\
    \midrule
    0.33 & ${0.01}^{+0.21}_{--0.21}$ & ${132.3}^{+2.3}_{-3.8}$ & 52.8783\\
    0.50 & ${0.012}^{+0.36}_{-0.012}$ & ${131.33}^{+4.3}_{-0.39}$ & 53.6116\\
    1.00 & ${0.18}^{+0.25}_{-0.18}$ & ${130.5}^{+2.2}_{-1}$ & 52.3628\\
    1.50 & ${0.32}^{+0.16}_{-0.12}$ & ${129.61}^{+1}_{-0.75}$ & 51.0397\\
    2.00 & ${0.34}^{+0.61}_{-0.31}$ & ${128.4}^{+2.5}_{-1.4}$ & 51.5343\\
    \bottomrule
    \end{tabular}


Filament 3
----------

ERRORS EVERYWHERE MY GOODNESS.

Closer look at error rootfinding 
================================

Filament 1, mu=0.33
-------------------
    Best fit B0 = 127.482 muG; eta2 = 5.678 (12 f'n calls)
    
    Bounded B0 below, in grid (78 function calls)
    Bounded B0 above, past grid (137 function calls)
        Bracket pts: B0 = 441.707 muG; eta2 = 3239.913
                     B0 = 494.366 muG; eta2 = 8991.098
        Error value: B0 = 485.162 muG; eta2 = 7452.972
    Bounded eta2 below, in grid (74 calls)
(!) Bounded eta2 above, in grid, PROBLEM: (54 calls for whole rigmarole)
        Started at: B0 = 127.482 muG; eta2 = 82.864
            but, rminarc too small (!).  Fit fails to explore larger B0 values
            after 4 calls, fit stops.  Code thinks we passed error threshold.
            
            So, in fact, we move backwards...
                B0 = 127.482 muG; eta2 = 68.665
                ...             ; eta2 = 68.665
                B0 = 213.132 muG; eta2 = 68.665
            This appears to be below the error threshold.  Now we pull out
            the rootfinder.  When checking the other point, we start at:
                B0 = 213.132 muG; eta2 = 82.864
                ...
                B0 = 222.162 muG; eta2 = 82.864
            Now, the code identifies this as inconsistent behavior!
            
            Code makes one more attempt to fix things, and searches between:
                B0 = 204.485 muG; eta2 = 56.899
                B0 = 231.577 muG; eta2 = 100.000
            NOPE, FAILURE.  Brentq fails and my code gives up, moves on.

_solution_: Generate new tables, and also take initial guess for "other param"
(B0 or eta2) from table?!


Filament 1, mu=0.5
------------------
    Best fit B0 = 115.956 muG; eta2 = 3.849
    Bound B0 below, in grid
(!) Bound B0 above, past grid, PROBLEM: sporadic box/resolution errors.
        No problem working up to B0 = 696.473 muG; eta2 = 62499.564
        Then, first instance: 
            B0 = 795.590 muG; eta2 = 163307.996
                Resolution error! at [ 2.]
        Code jumps down and keeps trying to increment eta2, but keeps erroring:
            B0 = 795.590 muG; eta2 = 79680.745
            ...
            B0 = 795.590 muG; eta2 = 119611.975
            B0 = 795.590 muG; eta2 = 190763.060
                Resolution error! at [ 2.]
            ...
            B0 = 795.590 muG; eta2 = 146893.065
                Box error! at [ 2.]
            ...
            Max eta2 that passes: B0 = 795.590 muG; eta2 = 146434.006
            Min eta2 that errors: B0 = 795.590 muG; eta2 = 146434.216
            
        Code stops at:
            B0 = 795.590 muG; eta2 = 146434.216 (which gives "Box error!")
        and thinks it has found the error crossing.
        It pulls out the brentq rootfinder and settles upon
            B0 = 774.021 muG; eta2 = 153164.184
                Resolution error! at [ 2.]
        But, this is patently incorrect.
    Bound eta2 below, in grid
(!) Bound eta2 above, past grid, PROBLEM -- same behavior as for mu=0.33
        (yes, brentq failed)

_solution_: address this box/resolution error, or set more conservative limits
            on eta2/B0.


Filament 1, mu=1.0
------------------
    Best fit B0 = 102.494 muG; eta2 = 2.590 (12 calls)
    Bound B0 below, in grid (62 calls)
    Bound B0 above, in grid (had to move backwards) (62 calls, coincidence...)
    Bound eta2 below, in grid (had to move backwards) (48 calls)
    Bound eta2 below, in grid (had to move backwards) (67 calls)

Filament 1, mu=1.5
------------------
    Best fit B0 = 97.081 muG; eta2 = 2.457
    Error bounding went smoothly

Filament 1, mu=2.0
------------------
    All went smoothly (it seems)


Filament 2, mu=0.33
-------------------
    Best fit B0 = 132.288 muG; eta2 = 0.015
(!) Bound B0 below, in grid:
        B0 = 128.533 muG; eta2 = 0.017 passes threshold
        Did not move on grid, checking backwards
            Settles on B0 = 128.533 muG; eta2 = 0.018
            Moves up to B0 = 129.271 muG; eta2 = 0.020
        Inconsistent fit behavior, tries to search between
            B0 = 341.068 muG; eta2 = 0.020 (where did this come from ?!?!?!)
            and
            B0 = 129.814 muG; eta2 = 104.465
        ...
        brentq fails, of course.  So we've got a BUG in my limits.
    Bound B0 above (833 function calls!!!!!), starts from:
        B0 = 341.068 muG; eta2 = 0.015  (an absurdly high value?!)
        ...
        B0 = 133.563 muG; eta2 = 0.152
        Finally reaches threshold
        B0 = 134.615 muG; eta2 = 0.219
(!) Bound eta2 below (1076 function calls!!!!!)
        (Annealing grid in 0.0 direction... WHAT?!)
        Starts at: B0 = 132.288 muG; eta2 = 0.000 and
        settles @: B0 = 132.818 muG; eta2 = 0.000
        "Hit grid edge 100.0"  CRAP.
        B0 = 132.818 muG; eta2 = 100.000.
        Now it pulls out `one_dir_root`, stepping all the way down to 0.
        Taking, forever.  But it eventually settles on a reasonable number.
        B0 = 134.575 muG; eta2 = 0.221
        THIS IS LARGER THAN OUR ORIGINAL eta2...
    Bound eta2 above, started at:
        B0 = 132.288 muG; eta2 = 100.000
            Keeps getting rminarc errors!  Too small...
            Eventually, steps backwards and finds a decent bound.

_solution_: code bug.  Cannot just catch IndexError because -1 is a valid index.
Addressed sgn = 0 edge case.


Filament 2, mu=0.5
------------------
(!) Best fit B0 = 131.330 muG; eta2 = 0.012  (suspicious; very close to starting
        values and fitting hit some intensity errors...)
    Bound B0 below, okay (B0 = 130.936 muG; eta2 = 0.044)
    Bound B0 above, okay (B0 = 135.655 muG; eta2 = 0.372)
(!) Bound eta2 below, hit limit (eta2 = 0.)
    Bound eta2 above, okay (B0 = 135.594 muG; eta2 = 0.376)

Filament 2, mu=1.0
------------------
    Best fit B0 = 130.545 muG; eta2 = 0.176
    Bound B0 below
        checking backwards... takes a while but root finder converges nicely
    Bound B0 above, okay and quick
(!) Bound eta2 below, hit limit eta2 = 0.
    Bound eta2 above, okay

Filament 2, mu=1.5
------------------
    Best fit B0 = 129.609 muG; eta2 = 0.316
(!) Bound B0 below -- finds inconsistent fit behavior
        Grid: B0 = 126.765 muG; eta2 = 0.088
        Move backwards: B0 = 126.765 muG; eta2 = 0.107
        Then, B0 = 127.385 muG; eta2 = 0.088
        Brentq now complains
            (I'm not sure what happened -- seems like it didn't step back.)
        Tries again w/ B0 = 126.744 muG; eta2 = 0.088
                   and B0 = 128.886 muG; eta2 = 0.088
        Eventually settles upon B0 = 128.861 muG; eta2 = 0.157
    Bound B0 above, looks okay
    Bound eta2 below, looks okay (though it takes a while, 286 function calls)
    Bound eta2 above, looks okay

_solution_: stupid bug. missing condition in while loop, so it didn't check
backwards correctly.  But, it got lucky here and found an okay value on second
attempt.

Filament 2, mu=2.0
------------------
    Best fit B0 = 128.390 muG; eta2 = 0.343
    Bound B0 below looks okay
(!) Bound B0 above
        Brackets between:
            B0 = 130.691 muG; eta2 = 0.873
            B0 = 130.937 muG; eta2 = 0.927
        Rootfinder:
            B0 = 130.691 muG; eta2 = 0.873
            B0 = 130.937 muG; eta2 = 0.907
        Rootfinder complains, values don't bracket any more.  Tries again with:
            B0 = 130.000 muG; eta2 = 0.718
            B0 = 131.089 muG; eta2 = 0.962
        This works out, settles upon B0 = 130.936 muG; eta2 = 0.927
        
        BUT, earlier we saw that this value was (likely) BELOW threshold
            B0 = 130.937 muG; eta2 = 0.907 (vs. 0.927)
        so we are slightly underestimating the error, I think...
    Bound eta2 below (takes a while)
    Bound eta2 above, okay

_solution_: sigh.  it will never be perfectly consistent.
if unsure, I can always follow up on it manually...
here my code is doing what it's supposed to, after all.
