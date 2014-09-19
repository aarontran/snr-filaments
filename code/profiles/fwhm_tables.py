"""
Code to generate tables of
1. FWHM, m_E measurements
2. best fit data

Read from specified data folders

Aaron Tran
September 2014
"""

from __future__ import division

import cPickle as pickle
from glob import glob
import numpy as np

from latex_table import LatexTable
import regdict_utils as rdu

# ===========
# FWHM tables
# ===========

def generate_fwhm_tab(rroot, kevs, labels):
    """LaTeX table of FWHM and m_E values from region dicts w/ FWHM data

    Inputs:
        rroot (str): filestem of region dictionaries to read
        kevs (str): energy band positions to compute m_E (should match labels)
        labels (list): energy band labels, strings
        title (str): table title
    Output: LatexTable object w/form:
        region, 0.7-1kev, 1-2kev, 2-7kev, 0.7-1kev, 1-2kev, 2-7kev
        1       fwhm      fwhm    fwhm    m_E       m_E     m_E
        2       ...
    """
    ltab = LatexTable(['Region'] + labels + labels,
                      ['{}'] + [2]*len(labels) + [1]*len(labels),
                      'FWHMs (replace w/ filament headers!)',
                      prec=2)

    n_regs = rdu.get_nregs(rroot)
    data_all = np.empty((n_regs, len(kevs))) * float('nan')
    m_e_all = np.empty((n_regs, len(kevs))) * float('nan')

    for rdict, n in rdu.regdict_load(rroot):
        data = np.array([rdict[lab]['fwhm'] for lab in labels])
        errs = [rdict[lab]['errs'] for lab in labels]  # neg, pos
        eps = np.mean(errs, axis=1)  # average of errors for calculations/fits
        inds = (np.where(np.isfinite(data)))[0]  # Consistent w/ data_init

        # Compute m_E
        kevs = np.array(kevs)  # to be sure
        me = m_expt(kevs, data)
        meeps = m_expt_err(kevs, data, eps)

        # Build row for latex table
        ltr_fwhm = ['{:d}'.format(n)]
        ltr_me = []

        # NOTE HUGE ASSUMPTION -- only first n indices are dropped
        # (i.e., once you find a FWHM in some band, you find a FWHM in all
        # bands above...)  Very fragile code.  Works for Tycho...
        idx = 0
        for m in xrange(len(labels)):
            if m not in inds:
                ltr_fwhm.extend([float('nan'),float('nan'),float('nan')])
                ltr_me.extend([float('nan'),float('nan')])
            else:
                ltr_fwhm.append(data[idx])
                ltr_fwhm.extend(errs[idx])
                if idx > 0:
                    ltr_me.append(me[idx-1])
                    ltr_me.append(meeps[idx-1])
                else:
                    ltr_me.extend([float('nan'),float('nan')])
            idx += 1

        ltr_all = ltr_fwhm + ltr_me  # Join FWHMs, m_E on one row
        ltab.add_row(*ltr_all)

        # Tack data onto aggregate for averaging
        data_all[n-1,:] = data
        m_e_all[n-1, 1:] = me

    ltr_all = ['avg']
    fwhms_avg = np.nanmean(data_all, axis=0)
    fwhms_std = np.nanstd(data_all, ddof=1, axis=0)
    m_e_avg = np.nanmean(m_e_all, axis=0)
    m_e_std = np.nanstd(m_e_all, ddof=1, axis=0)
    for fwhm, std in zip(fwhms_avg, fwhms_std):
        ltr_all.extend([fwhm, std, std])
    for m_e, std in zip(m_e_avg, m_e_std):
        ltr_all.extend([m_e, std])
    ltab.add_row(*ltr_all)

    return ltab

def m_expt(ens, widths):
    """Calculated following Ressler et al. [2014]
    Output array size is one less than that of ens/widths
    """
    e2, e1 = ens[1:], ens[:-1]
    w2, w1 = widths[1:], widths[:-1]
    return np.log(w2/w1) / np.log(e2/e1)

def m_expt_err(ens, widths, eps):
    """Propagate filament width errors in quadrature
    Output array size is one less than that of ens/widths
    """
    e2, e1 = ens[1:], ens[:-1]
    w2, w1 = widths[1:], widths[:-1]
    w2e, w1e = eps[1:], eps[:-1]

    me = m_expt(ens, widths)
    return np.abs(me * 1./np.log(e2 / e1) * np.sqrt((w1e/w1)**2 + (w2e/w2)**2))


if __name__ == '__main__':
    pass

