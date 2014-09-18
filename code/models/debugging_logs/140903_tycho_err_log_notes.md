Debugging, Tycho error calculations (Sept. 3)
=============================================

Engine 0: regions 1, 2, 3, 4
Engine 1: regions 5, 6, 7
Engine 2: regions 8, 9, 10
Engine 3: regions 11 (DEATH), (12), (13)

How did the third engine (2) finish so early, and where is my log output?!

1st engine log: ~16130 lines logged
2nd engine log: ~16600 lines logged
3rd engine log: ~8000 lines logged
4th engine log: ~3200 lines logged


In iPython notebook, print out the stdout log (the regular `display_output`).
Parse through with the regex expression:

    ^[ \t]+.*\n^[ \t]+.*\n^[^ \t\n]+.*$

In vim, escape the `+` quantifier:

    ^[ \t]\+.*\n^[ \t]\+.*\n^[^ \t\n]\+.*$

(this finds: 2 indented lines followed by non-indented line)


Region 11, mu = 1
-----------------

Okay, looking only at region 11 right now.
For mu = 1, best fit values are B0 = 2453.473, eta2=16118.971 ?!

Compare:
* mu = 0.00; B0 = 669.687, eta2 = 23.242
* mu = 0.33; B0 = 808.257, eta2 = 71.552
* mu = 0.50; B0 = 898.772, eta2 = 130.379
* mu = 1.00; B0 = 2453.473, eta2 = 16118.971 (?!?!?!)
* mu = 1.50; B0 = 416.951, eta2 = 6.906
* mu = 2.00; B0 = 378.196, eta2 = 4.929

The mu=1 fit started from B0 = 1300 microGauss and just blew up from there...
But, looking at the error finding -- it seems like that may genuinely be the
best fit.  How did error finding fail?

?! I'm surprised we didn't run into trouble w/ `one_dir_root_grid`.

1. Code anneals in negative B0 direction, starting at B0 = 310.822 muG
   but because I don't set eta2 intelligently, it's at 16118.971 still.
   BOX ERROR ENSUES.  Fit doesn't move anywhere, code thinks we found the error
   crossing.
   It starts moving backwards over 89 (!) grid points, starting from 0.
   Hitting box errors all the way, thinking the errors are too large.
   Until we hit B0 = 583.806 and the box error lifts for highest energy band
   fit is able to pull down eta2, and things look better.

   But it still can't find a crossing, up to B0=1352.978, eta2=1153.214
   I.e., even at B0=1352.978, it is able to compute all FWHMs.
   But, the chi-sqr is still too large.
   
   It moves the threshold back by one:
   idx = 109, idx\_adj = 110    (last valid index is 109)
   NOW, while loop catches that we hit the edge.
   check\_f\_adj = f(xgrid[109]) > 0 still, so it correctly reports we failed.

   THEN, at brentq limits, it tries to take:
   `x_in = xgrid[idx_adj]`
   `x_out = xgrid[idx]`
   and the code fails horribly.

In my new tables I "only" tabulated eta2 up to 10,000.

Solution
--------

See notes from 2014 September 9.
I basically throw out the error annealing on grid.
Now, just start from grid and let `one_dir_root` search for me.
Avoids all the ugly edge cases, and should be faster/simpler.


