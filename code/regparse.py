"""
Library functions for parsing DS9, CIAO region files
and assorted useful file IO (common to multiple scripts)

Aaron Tran
2014 June 23
(last modified: 2014 July 8)
"""

import os
import re

import ds9


def main():
    pass


######################
# Handy path functions
######################

def count_files_regexp(stem):
    """Number of grouped files with given directory stem, where
    the basename of the stem is a regex (use raw string as needed)
    """
    pattern = os.path.basename(stem)
    if os.path.dirname(stem) == '':
        files = [f for f in os.listdir('.') if re.match(pattern,f)]
    else:
        files = [f for f in os.listdir(os.path.dirname(stem))
                 if re.match(pattern, f)]
    return len(files)


def check_dir(stem, verbose=False):
    """Check if stem directory exists and create if needed"""
    stemdir = os.path.dirname(stem)
    if not os.path.isdir(stemdir) and stemdir is not '':
        if verbose:
            print 'stem {} in nonexistent directory'.format(stem)
            print 'Creating directory {}'.format(stemdir)
        os.makedirs(stemdir)


#################
# Region file I/O
#################

def load_ds9reg(fname, aux=False):
    """List of region strings from ds9 region file;
    headers checked for correct filetype (does not check coordsys)
    Trailing newlines are preserved
    """
    with open(fname, 'r') as f:
        header = f.readline()
        settings = f.readline()
        coordsys = f.readline()
        if 'DS9 version 4.1' not in header:
            print 'WARNING: potentially invalid region file'
            print 'First line was: {}'.format(header)
        regstrs = list(f)

    if aux:
        return regstrs, (header, settings, coordsys)
    return regstrs


def write_ds9reg(fname, regstrs, headers):
    """Inverse operation for load_ds9reg(..., aux=True)
    Does not check if file exists.
    """
    with open(fname, 'w') as f:
        for h in headers:
            f.write(h)
        for r in regstrs:
            f.write(r)


def conv_fk5_to_phys(f_in, f_img, f_out):
    """f_in (str) input ds9 region file, fk5 coords
    f_img (str) input FITS file on which regions are defined
    f_out (str) output ds9 region file, phys coords

    Warning: does NOT check for overwriting!
    """
    d = ds9.ds9()
    d.set('file ' + f_img)
    d.set('regions load ' + f_in)
    d.set('regions system physical')
    d.set('regions save ' + f_out)
    d.set('frame clear')  # In case of un-deletable regions
    d.set('exit')
    reload(ds9)  # Ad hoc fix



################
# Region parsing
################

def regparse(regstr):
    """Break a region specification into pieces.
    Returns three items: rtype (str), rnums (list of floats), rprops (str)
    Preserves trailing newlines"""
    rtype, rnums, rprops = re.split('[\(\)]', regstr)  # This assumes there are no other parentheses floating around...
    rnums = [float(x) for x in re.split('\s*,\s*', rnums)]   # Could improve by only matching first set of parentheses
    return rtype, rnums, rprops


if __name__ == '__main__':
    main()

