"""
Attempt at a test suite for simple/full model fitting

Run at terminal with $ nosetests models_all_exec_test.py -v

Aaron Tran
August 2014
"""

from __future__ import division

import cPickle as pickle
import numpy as np

from nose.tools import eq_, ok_
from nose.tools import assert_almost_equal as aeq_
from unittest import TestCase

import models_exec as mex
import snr_catalog as snrcat

# ==================
# SN 1006, test data
# ==================

KEVS = np.array([0.7, 1.0, 2.0]) 
DATA = np.array([35.5, 31.94, 25.34])
EPS = np.array([1.73, .97, 1.71])

# NOTE this absolute import is liable to break / be obsoleted, watch out
with open('tables/sn1006_grid-6-100-20_2014-07-25.pkl', 'r') as fpkl:
    TAB = pickle.load(fpkl)

# Table 7 (Filament 1), Ressler et al., 2014
TAB_7 = {0.0: {'eta2': (7.0, 4.0), 'B0': (165, 21)},
         1/3: {'eta2': (2.4, 0.9), 'B0': (130,  8)},
         0.5: {'eta2': (1.8, 0.6), 'B0': (123,  7)},
         1.0: {'eta2': (1.1, 0.4), 'B0': (113,  4)},
         1.5: {'eta2': (0.8, 0.3), 'B0': (108,  3)},
         2.0: {'eta2': (0.7, 0.3), 'B0': (105,  3)}}

# Table 8 (Filament 1), Ressler et al., 2014
TAB_8 = {0.0: {'eta2': (7.5, 2.0), 'B0': (142, 5)},
         1/3: {'eta2': (4.0, 1.3), 'B0': (120, 5)},
         0.5: {'eta2': (3.0, 1.1), 'B0': (112, 4)},
         1.0: {'eta2': (2.0, 1.0), 'B0': (100, 3)},
         1.5: {'eta2': (1.9, 1.2), 'B0': ( 95, 3)},
         2.0: {'eta2': (2.0, 1.0), 'B0': ( 92, 4)}}


# ==========
# Test cases
# ==========

class TestFullModelErrorsSN1006(TestCase):
    """Test Fitter() full model fits to SN 1006, Flmt 1"""

    @classmethod
    def setUpClass(self):
        print 'SETUP FITTER OBJ'

    @classmethod
    def tearDownClass(self):
        print 'TEARDOWN FITTER OBJ'



class TestSimpModelErrorsSN1006(object):
    """Test Fitter() simple model fits to SN 1006, Flmt 1"""

    @classmethod
    def setUpClass(self):
        snr = snrcat.make_SN1006()
        self.fobj = mex.Fitter(snr, KEVS, DATA, EPS, TAB)


    def test_fit_params(self):
        """Generate tests of best-fits w/ Sean's Table 7"""

        def check_pars(pars1, pars2, deltas, frac):
            """Check parameter equality within frac + within error"""
            for p1, p2, d in zip(pars1, pars2, deltas):
                aeq_(p1, p2, delta = min(frac*p1, frac*p2, d))

        for mu in TAB_7:
            dstr = ('Simple fit params agree w/ Sean to min(10%,error); '
                    'SN 1006, mu={:0.2f}'.format(mu))

            res = self.fobj.fitter_simp(mu)
            B0 = res.params['B0'].value
            eta2 = res.params['eta2'].value

            B0_tab, B0_err = TAB_7[mu]['B0']
            eta2_tab, eta2_err = TAB_7[mu]['eta2']

            check_pars.description = dstr
            p_1 = [B0*1e6, eta2]
            p_2 = [B0_tab, eta2_tab]
            err = [B0_err, eta2_err]

            yield check_pars, p_1, p_2, err, 0.1


    def test_fit_errors(self):
        """Simple fit: check best-fit errors correspond to d(chisqr) of +1.
        Rather poor design...
        
        need res, pstr to do everything
        need ci_val, dchisqr, chisqr_best for assertions"""

        def get_var_chisqr(pstr, x, res):
            res.params[pstr].value = x
            res.params[pstr].vary = False
            res.prepare_fit()
            res.leastsq()
            res.params[pstr].vary = True
            return res.chisqr

        def f_assert(x1, x2, x_true):
            """Small, dummy helper method for generator..."""
            aeq_(x1, x_true)
            aeq_(x2, x_true)

        for mu in TAB_7:
            
            res = self.fobj.fitter_simp(mu)
            chisqr_best = res.chisqr
            ci_dict = self.fobj.get_errs(res, 'simp', method='manual',
                                         ci_vals=(0.683,))

            for pstr in ci_dict:
                dstr = ('Simple fit errors have correct chisqrs; '
                        'SN 1006, mu={:0.2f}, {} errors'.format(mu, pstr))

                x1, x2 = np.sort([t[1] for t in ci_dict[pstr] if t[0]==0.683])
                cs1 = get_var_chisqr(pstr, x1, res)
                cs2 = get_var_chisqr(pstr, x1, res)
                f_assert.description = dstr
                yield f_assert, cs1, cs2, chisqr_best + 1.0

        # Manually resetting description
        f_assert.description = ('Simple fit errors have correct chisqrs; '
                                'SN 1006')


class TestErrorRootFinderGrid(TestCase):
    """Test mex.one_dir_root_grid"""
    def f_poly(self, x):
        return x**2 - 9

    @classmethod
    def setUpClass(self):
        xm = np.linspace(-10,10,100)
        xm2 = xm[np.abs(xm) < 2.7]  # all(f_poly(xm2) < 0), no crossings
        self.xmesh = xm
        self.xmesh_sm = xm2

    def test_crossings(self):
        """Error root finder on grid: +/- crossings on polynomial"""
        ind = np.searchsorted(self.xmesh, 0.5)

        pos_x = mex.one_dir_root_grid(self.f_poly, ind, +1, self.xmesh, True)
        neg_x = mex.one_dir_root_grid(self.f_poly, ind, -1, self.xmesh, True)

        aeq_(pos_x, np.searchsorted(self.xmesh, 3))  # Idx of 1st x > 3
        aeq_(neg_x, np.searchsorted(self.xmesh, -3) - 1)  # Idx of last x < -3

    def test_limits(self):
        """Error root finder on grid: crossings not found on grid"""
        ind = np.searchsorted(self.xmesh_sm, 0.5)

        pos_x = mex.one_dir_root_grid(self.f_poly, ind, +1, self.xmesh_sm, True)
        neg_x = mex.one_dir_root_grid(self.f_poly, ind, -1, self.xmesh_sm, True)

        eq_(pos_x, None)
        eq_(neg_x, None)

    def test_edge(self):
        """Error root finder on grid: start on grid edge"""
        pos_x = mex.one_dir_root_grid(self.f_poly, len(self.xmesh_sm)-1, +1,
                                      self.xmesh_sm, True)
        neg_x = mex.one_dir_root_grid(self.f_poly, 0, -1, self.xmesh_sm, True)

        eq_(pos_x, None)
        eq_(neg_x, None)

class TestErrorRootfinder(TestCase):
    """Test mex.one_dir_root"""

    def f_poly(self, x):
        return x**2 - 9

    def test_crossings(self):
        """Error root finder: +/- crossings on polynomial"""
        aeq_(3.0, mex.one_dir_root(self.f_poly, 0., 10., verbose=True))
        aeq_(-3.0, mex.one_dir_root(self.f_poly, 2.5, -10., verbose=True))

    def test_limits(self):
        """Error root finder: crossings not found"""
        eq_(1., mex.one_dir_root(self.f_poly, -2., 1., verbose=True))
        eq_(-2.9, mex.one_dir_root(self.f_poly, 2.9, -2.9, verbose=True))

    def test_start_on_crossing(self):
        """Error root finder: edge cases (search start on crossing, limit)"""
        # Initialize search ON crossing
        eq_(3., mex.one_dir_root(self.f_poly, 3., 4., verbose=True))
        # Initialize search AND limit on crossing
        eq_(3., mex.one_dir_root(self.f_poly, 3., 3., verbose=True))
        # Initialize search on limit (search direction ambiguous)
        eq_(10., mex.one_dir_root(self.f_poly, 10., 10., verbose=True))


if __name__ == '__main__':
    print 'Run me using nosetests!'

