README
2014 September 22

Fits to fake data generated with:

    fitter_fakeit = mex.Fitter(snrcat.make_tycho(), TYCHO_KEVS,
                               np.ones(len(TYCHO_KEVS)), np.ones(len(TYCHO_KEVS)),
                               None)  # data, eps, tab not needed

    # A range of parameters to try faking
    pars_fakeit = {1: {'mu': 1., 'eta2': 1., 'B0': 200e-6},
                   2: {'mu': 1., 'eta2': 1., 'B0': 500e-6},
                   3: {'mu': 2., 'eta2': 1., 'B0': 500e-6},
                   4: {'mu': 1., 'eta2': 0.01, 'B0': 200e-6},
                   5: {'mu': 0.5, 'eta2': 1., 'B0': 500e-6},
                   6: {'mu': 1., 'eta2': 100., 'B0': 500e-6}}

    TYCHO_DATA = {}
    for n, p in pars_fakeit.items():
        w = fitter_fakeit.width_full(p['mu'], p['eta2'], p['B0'])
        print w
        TYCHO_DATA[n] = TYCHO_KEVS, w, 0.05*w, (np.where(np.isfinite(w)))[0]

So, a small range of parameter space (varying mu, B0, eta2) and using FWHMs
with 5% errors.

    Full model: B0 = 200.000 muG; eta2 = 1.000; mu = 1.000; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 9.8186557   8.485788    6.56665655  5.75042416  5.11696599]
[ 8.17286016  7.05957576  5.4597385   4.78263781  4.2439443 ]
    Full model: B0 = 500.000 muG; eta2 = 1.000; mu = 1.000; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 2.53641162  2.20528697  1.67526515  1.45451828  1.29590852]
[ 2.05881785  1.77954854  1.37746807  1.20708466  1.07143607]
    Full model: B0 = 500.000 muG; eta2 = 1.000; mu = 2.000; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 2.37437855  2.03875823  1.53734967  1.35543134  1.22384801]
[ 1.91137499  1.63604656  1.25786951  1.11352188  1.01293559]
    Full model: B0 = 200.000 muG; eta2 = 0.010; mu = 1.000; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 9.22512153  7.6662725   5.38760228  4.41677362  3.63496033]
[ 7.65997282  6.364092    4.46760461  3.64440634  2.9866425 ]
    Full model: B0 = 500.000 muG; eta2 = 1.000; mu = 0.500; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 2.70305556  2.36141727  1.78500396  1.52947213  1.34087636]
[ 2.20550931  1.90988999  1.46875166  1.27258726  1.11049893]
    Full model: B0 = 500.000 muG; eta2 = 100.000; mu = 1.000; 
    init rminarc = [ 20.  20.  20.  20.  20.]; adapted = [ 7.27154171  6.85739795  6.15572641  5.79210124  5.45515902]
[ 6.05951743  5.71742325  5.12835301  4.82003129  4.53288337]
