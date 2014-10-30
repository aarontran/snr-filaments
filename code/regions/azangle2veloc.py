"""
Interpolate shock velocities given list of azimuth angles and measured
velocities

Used to interpolate shock velocities for Tycho's SNR

Aaron Tran
Fall 2014
"""

from __future__ import division

import argparse
from datetime import datetime
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(description=
             ('Interpolate shock velocities from azimuth angles; '
              'given discrete list of (az, veloc)'))

    parser.add_argument('azfile', help='Plaintext file of az angles')
    parser.add_argument('briantab', help=('Plaintext file of (az, v_shock) '
                                          'values; for Tycho, table from '
                                          'Williams et al. (2013)'))
    parser.add_argument('fout', help='Output filename (plaintext list of vs)')
    parser.add_argument('scale', help=('For Tycho, this is (3 kpc / 2.3 kpc) '
                                       '= 1.3043478260869565'))

    args = parser.parse_args()
    f_az, f_briantab = args.azfile, args.briantab
    f_out = args.fout
    scale = float(args.scale)

    az_angles = np.loadtxt(f_az)
    briantab = np.loadtxt(f_briantab)

    aztab = briantab[:,0]
    vstab = briantab[:,1] * 1e5  # Convert km/s to cm/s

    # Pad table on each end to catch edge cases while interpolating
    aztab = np.insert(aztab, [0, len(aztab)], [aztab[-1]-360., aztab[0]+360.])
    vstab = np.insert(vstab, [0, len(vstab)], [vstab[-1], vstab[0]])

    vs_out = np.empty(len(az_angles))

    # Interpolate linearly in Brian's table
    for n, az in zip(xrange(len(az_angles)), az_angles):
        az = az % 360
        ind = np.searchsorted(aztab, az)
        slope = (vstab[ind] - vstab[ind-1]) / (aztab[ind] - aztab[ind-1])
        vs_interp = vstab[ind-1] + slope * (az - aztab[ind-1])

        vs_out[n] = scale * vs_interp

    meta = 'Shock velocity values, scale factor {}\n'.format(scale)
    meta += 'Azimuth angles: {}\n'.format(os.path.abspath(f_az))
    meta += 'Table of (az, v_shock) values: {}\n'.format(os.path.abspath(f_briantab))
    meta += 'Generated: {}'.format(datetime.now())

    np.savetxt(f_out, vs_out, fmt='%.18f', header=meta)

    print 'Wrote interpolated shock velocities to {}'.format(f_out)


if __name__ == '__main__':
    main()
