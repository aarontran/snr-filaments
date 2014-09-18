"""
Script to prepare DS9 extracted profiles for FWHM fitting -- convert
x-coordinate to radial distance (arcsec), compute errors, compute fit_domain by
smoothing profiles and identifying downstream local minimum

Aaron Tran
July (June?) 2014, converted to .py from .ipynb Sept. 2014
"""

import argparse
import numpy as np

import regparse
import fsmooth

ACIS_PX2ARCSEC = 0.492  # Chandra ACIS pixel size, arcseconds
# http://cxc.harvard.edu/proposer/POG/html/chap6.html

def main():
    """Run from command line as, e.g.,
    python prepare_profile_fit.py 'regions-1.physreg' 'profiles/prf'
                                  'profiles/prf-cts' --bin 1
                                  --labels '0.7-1kev' '1-2kev' '2-7kev'
    Generates lots of files
    """
    parser = argparse.ArgumentParser(description=
             ('Make usable profile data for FWHM fitting.  Process output of '
              'ds9projplotter.py: convert to arcsec, compute errors. '
              'Outputs .npz and .dat files; downstream scripts use .npz, but '
              'plaintext (.dat) is nice for archiving.'))
    parser.add_argument('ds9physreg', help='DS9 region file, physical coords')
    parser.add_argument('inroot', help='Intensity profiles, file root')
    parser.add_argument('ctroot', help='Uncorrected count profiles, file root')
    parser.add_argument('outroot', help='Output profile root')

    parser.add_argument('-b', '--binsize', default=1, type=float,
                        help=('Bin size of images used to generate profiles; '
                        'set in CIAO\'s merge_obs, reproject_obs plists'))
    parser.add_argument('-l', '--labels', nargs='+',
                        help='Energy band labels for profiles; must be '
                             'sorted & consistent for subsequent analysis')

    parser.add_argument('-w', '--window', default='hanning',
                        help=('Window type for smoothing to set fit domain. '
                        'Options: flat, hanning, hamming, bartlett, blackman '
                        '(see fsmooth.std_smooth)'))
    parser.add_argument('-n', '--window-n', default=21, type=int,
                        help='Smoothing window length (pixels/points)')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose mode')

    args = parser.parse_args()
    ds9physreg = args.ds9physreg
    inroot, ctroot, outroot = args.inroot, args.ctroot, args.outroot
    binsize, labels = args.binsize, args.labels
    wtype, wlen = args.window, args.window_n
    verbose = args.verbose

    if not binsize:
        print 'No binsize supplied, assuming binsize 1 (0.492 arcsec/pixel)'
        print 'DEBUG: {}'.format(binsize)
        binsize = 1

    dims = get_proj_dims(ds9physreg)
    px2arcsec = binsize * ACIS_PX2ARCSEC

    fit_cuts = []

    for i in xrange(len(dims)):
        length, thickness = dims[i]  # Units: physical pixels
        thck = np.floor(thickness)  # Floor to be consistent w/ DS9,
                                    # for number of pixels integrated
        if verbose:
            print ('Region {:02d}: length={:0.2f}px, '
                   'thickness={:0.2f}px'.format(i+1, length, thck))

        reg_cuts = [i+1]  # First entry = region number

        for lab in labels:
            # Intensity data; x data is just 1, 2, 3, ...
            infile = '{0}_{1:02d}_band_{2}.dat'.format(inroot, i+1, lab)
            data = np.loadtxt(infile)
            x = (data[:,0] - 1) * px2arcsec  # Convert to arcsec, with x[0] = 0
            y = data[:,1]

            # Count data; counts already averaged over integration length
            ctfile = '{0}_{1:02d}_band_{2}.dat'.format(ctroot, i+1, lab)
            data_cts = np.loadtxt(ctfile)
            cts = data_cts[:,1]
            cts_err = np.sqrt(cts * thck) / thck  # Poisson: n/L +/- sqrt(n)/L
            cts_err[cts_err==0] = np.sqrt(1) / thck  # if n=0, err = sqrt(1)/L

            # Use count data to compute y errors
            np_errdict = np.seterr(divide='ignore', invalid='ignore')
            y_err = y * (cts_err / cts)
            # Edge cases if cts=0, or y=0 (cts_err addressed above)
            # if cts==0, use y +/- y (1 +/- 1 for Poisson error)
            # if y_err==0 (y==0), use minimum non-zero error
            y_err[cts==0] = y[cts==0]
            y_err[y_err==0] = np.amin(y_err[np.nonzero(y_err)])
            np.seterr(**np_errdict)  # Restore error warnings

            # Determine cuts in fit domain
            cut = fsmooth.ind_first_min(y,window_len=wlen,window=wtype)
            reg_cuts.append(cut)
            if verbose:
                tot_cts = np.sum(cts[cut:]) * thck
                avg_cts = np.mean(cts[cut:]) * thck
                print ('  Band {}: avg integrated cts {:0.3f}, total cts {} '
                       'in fit domain'.format(lab, avg_cts, tot_cts))

            # Save data to output files
            outfile_dat = '{0}_{1:02d}_band_{2}.dat'.format(outroot, i+1, lab)
            outfile_npz = '{0}_{1:02d}_band_{2}.npz'.format(outroot, i+1, lab)
            np.savetxt(outfile_dat, np.array((x, y, y_err)).T)
            np.savez(outfile_npz, x=x, y=y, y_err=y_err)

        fit_cuts.append(reg_cuts)

    fit_cuts = np.array(fit_cuts)

    # Save fit domain cuts to output files -- also note profile information
    cutfile_dat = '{0}_fit_cuts.dat'.format(outroot)
    cutfile_npz = '{0}_fit_cuts.npz'.format(outroot)
    hdr =  'Profile fit domain downstream cuts (array indices)\n'
    hdr += 'binsize={}, resolution={} arcsec\n'.format(binsize, px2arcsec)
    hdr += 'smoothing window: {}, length={}\n'.format(wtype, wlen)
    hdr += 'region number, ' + ', '.join(labels)
    np.savetxt(cutfile_dat, fit_cuts, fmt='%d', header=hdr)
    np.savez(cutfile_npz, cuts=fit_cuts)

    if verbose:
        print 'Wrote fit domain cuts to {}, {}'.format(cutfile_dat,cutfile_npz)
        print 'Done!'


def get_proj_dims(f_reg):
    """Get projection dimensions from region file, in physical coords

    Input: region filename to load
    Output: list of 2-tuples of length and thickness (radial dist, integration
            dist) in pixels
    """
    regstrs, headers = regparse.load_ds9reg(f_reg, aux=True)
    if 'physical' not in headers[2]:
        raise Exception('Expected physical coords, got {}'.format(headers[2]))

    dims = []
    for line in regstrs:
        rtype, rnums, rprops = regparse.regparse(line)
        if 'projection' not in rtype:
            raise Exception('Expected projection region; got {}'.format(rtype))
        x1,y1,x2,y2,w = rnums
        length = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        dims.append((length, w))
    return dims


if __name__ == '__main__':
    main()
