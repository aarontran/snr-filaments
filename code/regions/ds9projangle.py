"""
Compute azimuth angle between rim regions and remnant center

Aaron Tran
2014 August 6
"""

import argparse
from datetime import datetime
import numpy as np
import os

import regparse

RAD_2_DEG = 180. / np.pi


def main():
    """"""
    parser = argparse.ArgumentParser(description=('Get azimuth angles for rim '
                                     'regions relative to circle center.'))
    parser.add_argument('freg', help=('DS9 region file, physical coords, of '
                                      'projection regions around SNR'))
    parser.add_argument('fcirc', help=('DS9 region file, physical coords, of '
                                       'circle region centered on SNR'))
    parser.add_argument('fout', help='Output plaintext list of az angles')

    args = parser.parse_args()
    freg, fcirc, fout = args.freg, args.fcirc, args.fout

    p_circ = parse_circ(fcirc)

    az_vals = []
    rlist = regparse.load_ds9reg(freg)
    for rstr in rlist:
        p_proj = proj_center(rstr)
        az = RAD_2_DEG * get_az_angle(p_circ, p_proj)
        az_vals.append(az)

    # Save some information about how az angles were computed too
    meta = 'Source regions: {}\n'.format(os.path.abspath(freg))
    meta += 'Remnant circle: {}\n'.format(os.path.abspath(fcirc))
    meta += 'Angles generated: {}\n'.format(datetime.now())
    meta += 'Remnant center: ({:0.18f}, {:0.18f})'.format(*p_circ)

    np.savetxt(fout, az_vals, fmt='%.18f', header=meta)

    print 'Wrote azimuth angles (degrees CCW from N) to {}'.format(fout)


def get_az_angle(p1, p2):
    """Azimuth angle as measured from north, going east, with p1 at center.
    Range is [0, 2pi]"""
    x1, y1 = p1
    x2, y2 = p2
    xp, yp = x2-x1, y2-y1
    az = np.arctan2(-xp, yp)  # Measure CCW from north
    if az < 0:
        az = az + 2*np.pi
    return az


def proj_center(rstr):
    """Compute center coordinates of projection region.
    Partially lifted from ds9proj2ciao.py"""
    rtype, rnums, rprops = regparse.regparse(rstr)
    if 'projection' not in rtype:
        raise Exception('Invalid region, got {}'.format(rtype))
    x1, y1, x2, y2, t = rnums

    angle = np.arctan2(y2-y1, x2-x1)
    xm, ym = (x2+x1)/2, (y2+y1)/2  # Midpoint between p1, p2
    # Center point of box
    xc = xm - (t/2) * np.sin(angle)
    yc = ym + (t/2) * np.cos(angle)
    return xc, yc


def parse_circ(fcirc):
    """Get center coordinates of a circle from DS9 region file"""
    x0, y0 = -1, -1
    for rstr in regparse.load_ds9reg(fcirc):
        rtype, rnums, rprops = regparse.regparse(rstr)
        if 'circle' in rtype.lower() and x0 == -1:
            x0, y0 = rnums[:2]
    if x0 == -1:
        raise Exception('No circle found in circle-region file, liar.')
    return x0, y0


if __name__ == '__main__':
    main()
