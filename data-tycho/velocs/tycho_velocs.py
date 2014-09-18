import numpy as np
import matplotlib.pyplot as plt
from fplot import fplot

a = np.loadtxt('tycho_velocs.txt')

msk = a[:,1]>2500
nmsk = np.array(map(lambda x: not(x), msk))

plt.errorbar(a[:,0][msk], a[:,1][msk], yerr=a[:,2][msk], fmt='ob')
plt.errorbar(a[:,0][nmsk], a[:,1][nmsk], yerr=a[:,2][nmsk], fmt='or')
fplot('Azimuthal angle east of north (deg.)', r'Shock velocity (km s${}^{-1}$), with $D = 2.3$ kpc')
plt.show()

print 'Statistics of data with vs > 2500 km/s'
print 'Mean: {}, median: {}'.format(np.mean(a[:,1][msk]), np.median(a[:,1][msk]))
print 'Min: {}, max: {}'.format(min(a[:,1][msk]), max(a[:,1][msk]))