"""
Script to convert ds9 region files to CIAO region files
Aaron Tran
June 10, 2014
(last modified: June 19, 2014)

CIAO doesn't recognize boxes, so we must do the following:

1. read in DS9 file containing some projection regions
2. re-save DS9 file in physical coords
3. process physical coord DS9 file, and convert projection regions to boxes
4. load newly reprocessed DS9 file, then save as CIAO file
5. clean intermediate files (unless flagged otherwise)

Idea: this should be readily extended to help process any other
CIAO incompatible regions.
See translations here:
http://cxc.harvard.edu/ciao/ahelp/dmregions.html#DS9_regions
See script for converting ds9 regions with groups, to CIAO format:
http://cxc.harvard.edu/ciao/ahelp/dmgroupreg.html
"""

import argparse
import numpy as np
import os
import re

import ds9


def main():
    """Run from command line as python ds9proj2box.py -nv fnamein fnameout
    e.g.,
    python ds9proj2box.py -v ../data/2-7kev_mosaic.fits \
                             ../data/profiles_good.reg \
                             ../data/profiles_good_ciao.reg
    """
    
    parser = argparse.ArgumentParser(description=
             'Convert DS9 region files with projections to CIAO region files')
    
    parser.add_argument('img', help='FITS file on which regions are defined')
    parser.add_argument('input', help='Input DS9 region file (fk5 coords)')
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
        print 'Intermediate files:\n\t{0}\n\t{1}'.format(fname_phys, fname_physproc)
    
    # Convert ds9 file to physical coordinates
    d = ds9.ds9()  # I'm not sure why I'm getting the error
    # An instance of ds9 was found to be running before we could start the
    # 'xpans' name server.
    # It seems to work okay anyways.  Alternatively we can just stick to the
    # iPYthon notebook
    d.set('file ' + f_img)
    d.set('regions load ' + fname_in)
    d.set('regions system physical')
    d.set('regions save ' + fname_phys)
    d.set('frame clear')  # In case of un-deletable regions
    
    # Parse file in physical coordinates, write new ds9 region file
    with open(fname_phys, 'r') as f:
        # Check the file headers
        header = f.readline()
        settings = f.readline()
        coordsys = f.readline()
        if verbose:
            print 'Checking region file header...'
        if header != '# Region file format: DS9 version 4.1\n':
            print 'WARNING: potentially invalid region file!'
            print 'First line was: ' + header

        with open(fname_physproc, 'w') as fw:
            if verbose:
                print 'Writing regions to {}'.format(fname_physproc)
            fw.writelines([header, settings, coordsys])
            if verbose:
                print 'Parsing regions...'
            for line in f:
                if '# projection' in line:
                    fw.write(proj2box(line))
                else:
                    fw.write(line)
        if verbose:
            print 'Finished writing to {}'.format(fname_physproc)
        
    # Load processed region file
    d.set('file ../data/2-7kev_mosaic.fits')
    d.set('regions load ' + fname_physproc)
    # Save to CIAO format
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
    """Changes projection string (r) to box string.  No error checking"""
    rtype, rnums, rprops = re.split('[\(\)]', r)
    x1, y1, x2, y2, t = [float(x) for x in re.split(',', rnums)]

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
        return boxstr + ' #' + rprops  # projection doesn't match normal region formatting
    return boxstr + rprops


if __name__ == '__main__':
    main()
