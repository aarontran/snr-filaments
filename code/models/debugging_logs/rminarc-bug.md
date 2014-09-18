README on rminarc bug
=====================
2014 July 25

A bug in the full model code

Summary
-------
The f2py compiled Fortran code has a bug, for certain values of rminarc it
fails to compute the intensity profiles correctly (spitting out NaNs) and
returns invalid FWHMs.
This issue does not seem to occur in the plain Fortran code
(when compiled as `gfortran FullEfflength_mod.f` as usual).

The ad hoc solution I adopt for now is to guess and check rminarc,
and avoid changing rminarc in the tabulation code.

Quan M. Nguyen has suggested this might be a floating point error.

For reference, some software versions:

> `GNU Fortran (Homebrew gcc 4.8.3_1) 4.8.3`
> `Python 2.7.7 (default, Jun 18 2014, 09:51:23)`
> `np.__version__ == '1.8.1'`

Observations and notes
----------------------
1. certain numbers are susceptible, and they are seemingly random
   e.g. 73.71, but not 73.72
2. Error occurs when f2py module is called from tabulating functions,
   called from `models_all.width_cont` wrapping function,
   and when called from f2py compiled `fullmodel.so` directly.
   Hence my belief that the error occurs at f2py compile time.
3. It appears to be independent of the parameters B0, eta2, mu, although
   I have not fully vetted this.

Example: call from tabulating code.

    Function call with B0 = 93.871 muG; eta2 = 0.010; mu = 0.333
    500000000.0 125000000.0 2.96e+19 900 2.2 73.71 1 400 100 500
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at  0.69999999999999996     
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   1.0000000000000000     
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   2.0000000000000000     
        Model fwhms = [ 900.  900.  900.]
        
Call `width_cont` directly from command line

    p = lmfit.Parameters()
    p.add('B0', value=93.871e-6)
    p.add('eta2', value=0.01,vary=False)
    p.add('mu',value=1./3,vary=False)

    ma.width_cont(p, ma.SN1006_KEVS, ma.snrcat.make_SN1006(), rminarc=73.71)
        Function call with B0 = 93.871 muG; eta2 = 0.010; mu = 0.333
    [ 0.7  1.   2. ] 9.3871e-05 0.01 0.333333333333
    500000000.0 125000000.0 2.96e+19 900 2.2 73.71 1 400 100 500
     Rim falloff towards shock is weird?

Call from precompiled `fullmodel.so`:

    fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 500000000.0,125000000.0,2.96e+19,900,2.2,73.71,1,400,100,500)
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at  0.69999999999999996
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   1.0000000000000000
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   2.0000000000000000

Try changing the value of rminarc, holding all other parameters constant.
Works for 60, 73, 58.5, 73.5, 73.8, 73.72, 73.5.
Fails on 74, 73.71.  Some example outputs:

    In [52]: fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 5e8, 1.25e8,2.96e19,900,2.2,73.5,1,400,100,500)
    Out[52]: array([ 35.83,  29.95,  21.13])

    In [54]: fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 5e8, 1.25e8,2.96e19,900,2.2,73.71,1,400,100,500)
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at  0.69999999999999996     
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   1.0000000000000000     
     Rim falloff towards shock is weird?
     Box Length Error (xmin) at   2.0000000000000000     
    Out[54]: array([ 900.,  900.,  900.])

    In [55]: fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 5e8, 1.25e8,2.96e19,900,2.2,73.72,1,400,100,500)
    Out[55]: array([ 35.75,  29.86,  21.01])

The original fortran code doesn't appear to display this error.

     Enter B0
    9.3871e-05
     Enter Eta2
    0.01
     Enter mu
    0.3333333333
     Enter rminarc
    73.71
      0.69999999999999996        1.0000000000000000        2.0000000000000000        3.0000000000000000        9.3870999999999996E-005   1.0000000000000000E-002  0.33333333329999998     
       500000000.00000000        125000000.00000000        2.9600000000000000E+019   900.00000000000000        2.2000000000000002        73.709999999999994                1         400         100         500
      0.69999999999999996       keV:    35.749349999999971     
       1.0000000000000000       keV:    29.852549999999965     
       2.0000000000000000       keV:    21.007350000000013     
       3.0000000000000000       keV:    16.953299999999992

One more retry, with f2py modified to spit out numbers

    import fullmodel as fm
    fm.readfglists()
    fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 5e8, 1.25e8,2.96e19,900,2.2,73.71,1,400,100,500)

When I examine the output intensity distributions, the ones from f2py
are full of NaNs (for rminarc=73.71).  Electron distributions look okay.
And, again, it works just fine for rminarc=73.72.

    fm.fullefflengthsub([0.7,1.,2.], 9.3871e-05, 0.01, 0.3333333333, 5e8, 1.25e8,2.96e19,900,2.2,73.72,1,400,100,500)
      0.69999999999999996        1.0000000000000000        2.0000000000000000        9.3870999999999996E-005   1.0000000000000000E-002  0.33333333329999998     
       500000000.00000000        125000000.00000000        2.9600000000000000E+019   900.00000000000000        2.2000000000000002        73.719999999999999                1         400         100         500
    Out[7]: array([ 35.7542,  29.8566,  21.0102])

One more remark:
I ran the tabulation code on rminarc=73.71, but with different parameters:
eta2 = 0, B0 = 150, mu = 0.  Which leads me to suspect that the bug should be
parameter independent.  And, when I ran the code the first time, I didn't run
into this bug (not seen in printed output)
