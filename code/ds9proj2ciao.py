"""
Convert ds9 region files to CIAO format; manually convert projections to boxes

Aaron Tran
2014 June 10
(last modified: 2014 July 8)

1. read DS9 file containing projection regions
2. save DS9 file in physical coords.
3. read DS9 file with phys. coords. and convert projections to boxes,
   save DS9 file with boxes in phys coords.
4. convert DS9 file to CIAO file
5. clean intermediate files (unless flagged otherwise)

Idea: this should be readily extended to help process any other
CIAO incompatible regions.
Region translations:
http://cxc.harvard.edu/ciao/ahelp/dmregions.html#DS9_regions
Script for converting ds9 regions with groups to CIAO format:
http://cxc.harvard.edu/ciao/ahelp/dmgroupreg.html
"""

import argparse
import numpy as np
import os
import re

import ds9

import regparse


def main():
    """Run from command line as python ds9proj2ciao.py -v fnamein fnameout"""

    parser = argparse.ArgumentParser(description=
             'Convert DS9 region files with projections to CIAO region files')
    
    parser.add_argument('img', help='FITS file on which regions are defined')
    parser.add_argument('input', help='Input DS9 region file (any coords)')
    parser.add_argument('output', help='Output CIAO region file')
    parser.add_argument('-n', '--noclean', action='store_true',
                        help='Do not remove intermediate data files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    
    # Parse arguments
    args = parser.parse_args()
    f_img = args.img
    fname_in, fname_out = args.input, args.output
    noclean, verbose = args.noclean, args.verbose
    
    # Intermediate filenames
    fname_phys = fname_in + '.phys'
    fname_physproc = fname_phys + '.proc'
    if verbose:
        print 'Verbose mode enabled'
        print 'Intermediate files:\n\t{0}\n\t{1}'.format(fname_phys,
                                                         fname_physproc)
    
    # Convert ds9 file to physical coordinates
    regparse.conv_fk5_to_phys(fname_in, f_img, fname_phys)
    
    # Parse file in physical coordinates, write new ds9 region file
    regstr_phys, headers = regparse.load_ds9reg(fname_phys, aux=True)

    regstr_physproc = []
    if verbose:
        print 'Converting projections to boxes'
    for line in regstr_phys:
        if '# projection' in line:
            regstr_physproc.append(proj2box(line))
        else:
            regstr_physproc.append(line)

    if verbose:
        print 'Writing regions to intermediate file {}'.format(fname_physproc)
    regparse.write_ds9reg(fname_physproc, regstr_physproc, headers)

    # Load processed region file and save to CIAO format
    d = ds9.ds9()
    d.set('file ' + f_img)
    d.set('regions load ' + fname_physproc)
    d.set('regions system physical')
    d.set('regions format ciao')
    d.set('regions save ' + fname_out)
    if verbose:
        print 'CIAO output written to {}'.format(fname_out)    
    
    # Cleanup
    if not noclean:
        os.remove(fname_phys)
        os.remove(fname_physproc)
        if verbose:
            print 'Intermediate files removed'
    d.set('exit')
    if verbose:
        print 'Done!'


def proj2box(r):
    """Changes projection string to box string."""
    rtype, rnums, rprops = regparse.regparse(r)
    x1, y1, x2, y2, t = rnums

    if 'projection' not in rtype:
        raise Exception('Invalid region, got {}'.format(rtype))

    # Calculate various box parameters
    w = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    h = t
    angle = np.arctan2(y2-y1, x2-x1)
    xm, ym = (x2+x1)/2, (y2+y1)/2  # Midpoint between p1, p2
    # Center point of box
    xc = xm - (t/2) * np.sin(angle)
    yc = ym + (t/2) * np.cos(angle)

    boxstr = 'box({0},{1},{2},{3},{4})'.format(xc, yc, w, h,angle*180/np.pi)
    if rprops is not '\n':
        return boxstr + ' #' + rprops  # projection formatting is unusual
    return boxstr + rprops


if __name__ == '__main__':
    main()

