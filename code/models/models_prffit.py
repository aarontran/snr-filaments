"""
Short one-off script to try profile fitting for radio,X-ray w/ new code/data
Aaron Tran
November 2014
"""

import matplotlib.pyplot as plt
import numpy as np

import models
import snr_catalog as snrcat

# Because I'm just prototyping at this point
from IPython import get_ipython
ipython = get_ipython()
_wd = ipython.magic('%pwd')
ipython.magic('%cd ~/snr-research/VLA')
import plot_radio_xray_prfs as prxp
ipython.magic('%cd {}'.format(_wd))

def main():
    """Generate some fits.  Some stuff from interactive sessions"""
    prf_radio, prf_xhard = load_data(1)
    snr = snrcat.make_tycho()
    kevs_radio = np.array([5.6865e-9])
    kevs_xhard = np.array([4.5])  # Different than what we usually do

    # RADIO
    res_r = models.profile_fit(snr, kevs_radio, [prf_radio[:,1]],
            [np.ones_like(prf_radio[:,1])], [prf_radio[:,0]],
            1.0, eta2=1.0, B0=50e-6, r_trans=0, amp=1, ab_fit=0.005,
            mu_free=False, eta2_free=False, model_kws={'idamp':True})

    # XRAY
    res_x = models.profile_fit(snr, kevs_xhard, [prf_xhard[:,1]],
            [np.ones_like(prf_xhard[:,1])], [prf_xhard[:,0]],
            1.0, eta2=1.0, B0=25e-6, r_trans=0, amp=1, ab_fit=0.004,
            mu_free=False, eta2_free=False, model_kws={'idamp':True})

    # ATTEMPT TO PLOT THINGS
    plt.ion()
    plt_radio(res_r)
    plt_xhard(res_x)

# Data munging

def load_data(num):
    """Load preliminary data according to num
    num=1,2 (NNE,NW) implemented
    """
    prfs, scales = prxp.load_prf_data(num, pref='~/snr-research/VLA')
    prf_radio, _, prf_xhard = prfs

    if num == 1:
        msks_radio = (prf_radio[:,1] > 0,
                      prf_radio[:,0] > 230)
        msks_xhard = (prf_xhard[:,0] > 228,
                      prf_xhard[:,0] < 242)
    elif num == 2:
        msks_radio = (prf_radio[:,1] > 0,
                      prf_radio[:,0] > 232)

        msks_xhard = (prf_xhard[:,1] > 1e-9,
                      prf_xhard[:,0] > 228,
                      prf_xhard[:,0] < 240)

    msk_radio = reduce(np.logical_and, msks_radio)
    msk_xhard = reduce(np.logical_and, msks_xhard)

    return prf_radio[msk_radio], prf_xhard[msk_xhard]

# Plotting

def plt_radio(p):
    plt.gca().clear()
    _ = models.width_cont(p, kevs_radio, snr, rminarc=20,
                          idamp=True, damp_ab=p['ab'].value,
                          irad_adapt=False, get_prfs=True)
    _, igrids, rgrids = _

    plt.plot(rgrids[0], igrids[0], '-k.')
    plt.plot(prf_radio[:,0] + p['r_trans'].value,
             prf_radio[:,1]*p['amp'].value, 'or')
    plt.axvspan(240*(1-p['ab'].value), 240, alpha=0.1, color='green')
    plt.ylim(ymin=0)
    plt.show()

def plt_xhard(p):
    plt.gca().clear()
    _ = models.width_cont(p, kevs_xhard, snr, rminarc=20,
                          idamp=True, damp_ab=p['ab'].value,
                          irad_adapt=False, get_prfs=True)
    _, igrids, rgrids = _

    plt.plot(rgrids[0], igrids[0], '-k.')
    plt.plot(prf_xhard[:,0] + p['r_trans'].value,
             prf_xhard[:,1]*p['amp'].value, 'or')
    plt.axvspan(240*(1-p['ab'].value), 240, alpha=0.1, color='green')
    plt.ylim(ymin=0)
    plt.show()


if __name__ == '__main__':
    main()
