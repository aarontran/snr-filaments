"""

Interpolate either shock velocities or proper motions, given a list of azimuth
angles (generally, any azimuth-angle dependent parameter)

Shock velocities: intended for data of Williams+(2013), which in turn averages
data from Katsuda+(2010) and Reynoso+(1997) (X-ray, radio proper motion
respectively)

Proper motions: intended for data of Katsuda et al. (2010), using 2000-2007
baseline wherever possible (2003-2007 baseline in SSW of remnant, that was
obscured in 2000 observation)

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
             ('Linearly interpolate shock velocities or proper motions '
              'from azimuth angles.  Assumes one or the other.'
              'Outputs plaintext file w/ list of interpolated values'))

    parser.add_argument('azfile', help='Plaintext file of az angles')
    parser.add_argument('tabfile', help=('Plaintext table of az angles & qty '
                                         'to interpolate in 1st/2nd columns'))
    parser.add_argument('fout', help='Output plaintext filename')
    parser.add_argument('-b', '--briantab', action='store_true',
                        help='Flag if using Brian\'s velocs')
    parser.add_argument('-s', '--scale', default=3./2.3, type=float,
                        help=('Scale factor if using Brian\'s velocs; default '
                              'value is (3 kpc / 2.3 kpc) = 1.304348'))
    parser.add_argument('-n', '--nearest', action='store_true',
                        help='Use nearest nhbr instead of linear interp')

    args = parser.parse_args()
    f_az, f_tab, f_out = args.azfile, args.tabfile, args.fout
    briantab, scale = args.briantab, args.scale
    nearest = args.nearest

    az_angles = np.loadtxt(f_az)
    tab = np.loadtxt(f_tab)
    az_tab, par_tab = tab[:,0], tab[:,1]
    if briantab:
        par_tab *= 1e5  # Convert km/s to cm/s

    # Pad table on each end to catch edge cases while interpolating
    az_tab = np.insert(az_tab, [0, len(az_tab)],
                       [az_tab[-1]-360., az_tab[0]+360.])
    par_tab = np.insert(par_tab, [0, len(par_tab)],
                        [par_tab[-1], par_tab[0]])

    par_out = np.empty(len(az_angles))

    # Interpolate linearly in table
    # handles edge cases + unsorted az_angles
    for n, az in zip(xrange(len(az_angles)), az_angles):
        az = az % 360
        ind = np.searchsorted(az_tab, az)  # Should not be zero-th/last index

        if not nearest:  # Default, linear interpolation
            slope = (par_tab[ind] - par_tab[ind-1]) / (az_tab[ind] - az_tab[ind-1])
            par_out[n] = par_tab[ind-1] + slope * (az - az_tab[ind-1])
        else:  # Nearest neighbor interpolation
            az_left, az_right = az_tab[ind-1], az_tab[ind]
            ind_nearest = ind-1 if abs(az_left - az) < abs(az_right - az) else ind
            par_out[n] = par_tab[ind_nearest]

        if briantab:
            par_out[n] *= scale

    if briantab:
        meta = 'Shock velocity values, scale factor {}\n'.format(scale)
    else:
        meta = 'Proper motion values\n'
    meta += 'Azimuth angles: {}\n'.format(os.path.abspath(f_az))
    meta += 'Table of interp values: {}\n'.format(os.path.abspath(f_tab))
    if nearest:
        meta += 'Interpolation method: nearest\n'
    else:
        meta += 'Interpolation method: linear\n'
    meta += 'Generated: {}'.format(datetime.now())

    np.savetxt(f_out, par_out, fmt='%.18f', header=meta)

    print 'Wrote interpolated values to {}'.format(f_out)


if __name__ == '__main__':
    main()
