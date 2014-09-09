"""
Code to work with full model code interactively

Aaron Tran
September 2014
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

import lmfit

def main():
    pass

def manual_err(fobj, conf_intv):
    """Interactively find errors for (fobj or user-provided) best fit

    fobj.get_fit(...)
    for each free parameter:
        (i) user fixes one param's value
        (ii) user inputs guesses for all remaining free params
        (iii) compute chi-sqr value (whether from user guesses, or a fit
             starting from user guesses)
        (iv) if chi-sqr value is LESS than error threshold,
              go back to (i) and choose new fixed value
        (v) if chi-sqr value is GREATER Than error threshold,
            give user choice to 1. give new guesses for free params,
            2. try refitting, 3. go back to (1) and choose new fixed value,
        (vi) in BOTH cases, give user option to STOP and keep this best
             param
        (vii) also do this TWICE to get upper/lower bounds
    return error bounds...
    """


def manual_fit(snr, kevs, data, eps):
    """Interactively prompts user for values of B0, eta2, mu
    Prints FWHMs, residuals, and chi-squared values in the process
    And, plots the data + fitted model
    Also prints m_E, computed point to point (just for comparison)

    Note: doesn't really deal with bad input

    Input: snr object and data (x, y, eps) to be fitted
    Output: lmfit.Parameters object with best obtained fit values
    """

    print 'Manual fitting routine: enter q to print best fit and quit.\n'

    chisq_best = float('inf')
    w_model_best = None
    p_best = None

    while True:
        B0 = _get_float('Enter B0 (G): ')
        eta2 = _get_float('Enter eta2 (-): ')
        mu = _get_float('Enter mu (-): ')

        if any(np.isnan([B0, eta2, mu])):
            break  # Received exit/quit request

        p = lmfit.Parameters()
        p.add('B0', value=B0)
        p.add('eta2', value=eta2)
        p.add('mu', value=mu)

        adj = _get_float('Change model settings? Enter 0/1 for no/yes: ')
        if adj == 0:
            w_model = width_cont(p, kevs, snr)
        else:
            rminarc = _get_float('rminarc (default: 60): ')
            icut = _get_float('icut (default: 1): ')
            # Could add option to change resolution too
            w_model = width_cont(p, kevs, snr, rminarc=rminarc, icut=icut)

        chisq = chi_squared(data, eps, w_model)
        w_resid = (w_model - data)/eps

        print ('\nObsvd fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*data)
        print ('Model fwhms: ' + len(kevs)*'{:0.2f}, ')[:-2].format(*w_model)
        print ('Wghtd resid: ' + len(kevs)*'{:0.3f}, ')[:-2].format(*w_resid)
        print 'Chi^2: {:0.3f}'.format(chisq)

        plt.clf()
        plt.errorbar(kevs, data, eps, fmt='bo')
        plt.plot(kevs, w_model, '-k.')
        plt.draw()
        plt.show(block=False)

        if chisq < chisq_best:
            print '\nImproved fit, saving parameters'
            chisq_best = chisq
            w_model_best = w_model
            p_best = p
        else:
            print '\nBest fit so far:'
            print ('Best fwhms:  ' +
                   len(kevs)*'{:0.2f}, ')[:-2].format(*w_model_best)
            print 'Chi^2: {:0.3f}'.format(chisq_best)

        print '--------------------------------'
    
    if p_best is not None:
        plt.close()
        print '\nDone with manual fit. Best fit parameters are:'
        for key in p_best.keys():
            print '{} = {:g}'.format(p_best[key].name, p_best[key].value)

    return p_best


def _get_float(prompt):
    """Prompt user to input float, checking for exits / bad floats"""
    while True:
        try:
            uinput = raw_input(prompt).strip().lower()
            if uinput in ['q', 'quit', 'exit']:
                uinput = float('NaN') # Pass NaN to indicate exit/quit
            else:
                uinput = float(uinput)
            break  # Exit loop and return
        except ValueError:
            print '\nInvalid float, try again (enter q to quit)\n'
    return uinput


if __name__ == '__main__':
    main()
