"""
One-off script to plot proper motions from Katsuda et al. (2010)
Aaron Tran
"""

import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('tycho-motion.txt')

def bin_f(vec, a, b, f=np.mean):
    """For vector with columns of x, y coords, get mean y for a<x<b"""
    msk = np.logical_and(vec[:,0] >= a, vec[:,0] < b)
    if any(msk):
        return f(vec[msk,1])
    return np.nan

bin_left = np.arange(0, 360, 60)
bin_right = np.arange(60, 420, 60)
bin_ctr = np.arange(30, 360, 60)
a_mean = map(lambda x1,x2: bin_f(a,x1,x2,f=np.mean), bin_left, bin_right)
a_std = map(lambda x1,x2: bin_f(a,x1,x2,f=np.std), bin_left, bin_right)

plt.errorbar(a[:,0], a[:,1], yerr=a[:,2], fmt='ob')
plt.errorbar(bin_ctr, a_mean, xerr=30*np.ones(len(a_mean)),
             yerr=a_std, fmt='or')

plt.xlabel('Azimuth angle east of north (deg.)')
plt.ylabel('FS proper motion (arcsec/yr)')
plt.tight_layout()
plt.xlim((0,360))
plt.savefig('tycho-motions.pdf')
plt.show()

print 'Data statistics'
print 'Mean: {}, median: {}'.format(np.mean(a[:,1]), np.median(a[:,1]))
print 'Min: {}, max: {}'.format(min(a[:,1]), max(a[:,1]))

hdr = ('Binned (60 deg. sectors) Tycho proper motions; Katsuda et al. (2010)\n'
       'Column 1: bin centers (az angle in deg., east of north)\n'
       'Column 2: proper motion (arcsec/yr)\n'
       'Column 3: stdev of averaged proper motion data (arcsec/yr)')
np.savetxt('tycho-motion-bin.txt', np.array([list(bin_ctr), a_mean, a_std]).T,
           fmt='%.4f', header=hdr)
