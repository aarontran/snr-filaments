"""
Plot roughly obtained X-ray/radio profiles

# Note: scale/cmap settings not right for 1-1.7 keV image

ds9 -rgb \
    -red ../VLA/TYCHO_IF1.FITS \
        -scale limits 0.00015 0.0035 \
        -cmap value 2.25 0.55 \
    -green 1-1.7kev_mosaic.fits \
        -scale limits 7e-9 4e-7 \
        -cmap value 3.59 0.12 \
    -blue 4-7kev_mosaic.fits \
        -scale limits 7e-9 4e-7 \
        -cmap value 3.59 0.12 \
    -rgb lock scale yes \
    -asinh

    pan: 0:25:15.764, +64:08:25.00

"""
import numpy as np
import matplotlib.pyplot as plt

nyrs = 2009 - 1994  # Time between Chandra and VLA observations
arcsec_per_yr = 0.3  # Katsuda et al. 2010
px2arcsec_radio = 0.5  # Image scaling
px2arcsec_xray = 0.492
r_shock = 240  # arcsec
r_rise = 2  # Dist of x-ray peak behind shock (arbitrary guess!)

r_min = 200  # arcsec, plot bounds
r_max = 250

nums = [1,2,3,4]  # Labeling
names = ['NNE', 'NW', 'WNW', 'SW']


def main():
    """Make subplot grid of profiles for all regions"""
    plt.figure(figsize=(8,6))

    for num, name in zip(nums, names):

        prfs, scales = load_prf_data(num)
        prf_radio, prf_xsoft, prf_xhard = prfs
        scale_radio, scale_xsoft, scale_xhard = scales

        plt.subplot(2, 2, num)
        plt.plot(prf_radio[:,0], prf_radio[:,1] / scale_radio, '-r.',
                 label='1.375 GHz')
        plt.plot(prf_xsoft[:,0], prf_xsoft[:,1] / scale_xsoft, '-g.',
                 label='1-1.7 keV')
        plt.plot(prf_xhard[:,0], prf_xhard[:,1] / scale_xhard, '-b.',
                 label='4-7 keV')

        plt.xlim(r_min,r_max)
        plt.ylim(ymin=0)
        if num == 3:
            plt.ylim(ymax=1.5)

        # Annotations
        if num == 1:
            plt.legend(loc='best')
        if num in (3,4):
            plt.xlabel('Radial distance (arcsec)')
        if num in (1,3):
            plt.ylabel('Intensity (arbitrary scalings)')

        xpos_name = 0.58 if num == 1 else 0.95
        plt.text(xpos_name, 0.87, name, transform=plt.gca().transAxes,
                 horizontalalignment='right')

    plt.savefig('prf-flmts-subplot.pdf')
    plt.show()


def load_prf_data(num, pref=None):
    """Load data from 3 bands, shift and scale"""
    if pref is None:
        prf_radio = np.loadtxt('prf-flmts/prf_{:02d}_band_1.375GHz.dat'.format(num))
        prf_xsoft = np.loadtxt('prf-flmts/prf_{:02d}_band_1-1.7keV.dat'.format(num))
        prf_xhard = np.loadtxt('prf-flmts/prf_{:02d}_band_4-7keV.dat'.format(num))
    else:
        prf_radio = np.loadtxt(pref+'/prf_{:02d}_band_1.375GHz.dat'.format(num))
        prf_xsoft = np.loadtxt(pref+'/prf_{:02d}_band_1-1.7keV.dat'.format(num))
        prf_xhard = np.loadtxt(pref+'/prf_{:02d}_band_4-7keV.dat'.format(num))

    prfs_all = (prf_radio, prf_xsoft, prf_xhard)

    # Translate and align data by radial coordinate

    prf_radio[:,0] *= px2arcsec_radio  # Convert to arcseconds
    prf_xsoft[:,0] *= px2arcsec_xray
    prf_xhard[:,0] *= px2arcsec_xray

    prf_radio[:,0] += nyrs * arcsec_per_yr  # Proper motion (arcsec)

    ind_pk_init = np.argmax(prf_xhard[:,1])  # 4-7 keV peak
    r_pk_init = prf_xhard[ind_pk_init, 0]

    def rescale_to_rs(prf):  # Modifies input
        prf[:,0] += (r_shock - r_pk_init - r_rise)
        return prf
    def mask_r_domain(prf):
        r_mask = np.logical_and(r_min < prf[:,0], prf[:,0] < r_max)
        return prf[r_mask]

    map(rescale_to_rs, prfs_all)
    prf_radio, prf_xsoft, prf_xhard = map(mask_r_domain, prfs_all)

    # Rescale intensity values
    #del ind_pk_init, r_pk_init  # must update for new coordinates
    ind_pk = np.argmax(prf_xhard[:,1])
    r_pk = prf_xhard[ind_pk, 0]

    ind_pk_radio = np.searchsorted(prf_radio[:,0], r_pk)
    ind_pk_xsoft = np.searchsorted(prf_xsoft[:,0], r_pk)

    scale_radio = np.amax(prf_radio[(ind_pk_radio-10):, 1])  # Search 10 pts behind \
    scale_xsoft = np.amax(prf_xsoft[(ind_pk_xsoft-10):, 1])  # the 4-7 keV peak
    scale_xhard = prf_xhard[ind_pk, 1]
    scales_all = (scale_radio, scale_xsoft, scale_xhard)

    return prfs_all, scales_all


if __name__ == '__main__':
    main()
