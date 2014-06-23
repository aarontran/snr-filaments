"""
Library functions for parsing DS9, CIAO region files
and assorted useful file IO (common to multiple scripts)

Aaron Tran
2014 June 23
"""

import os
import re


def main():
    print 'Nothing to see here...'


#################
# Region file I/O
#################

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


def load_ds9reg(fname):
    """Returns list of region strings for ds9 region, with headers checked for correctness"""
    with open(fname, 'r') as f:
        header = f.readline()
        settings = f.readline()
        coordsys = f.readline()
        if 'DS9 version 4.1' not in header:
            print 'WARNING: potentially invalid region file'
            print 'First line was: {}'.format(header)
        if 'physical' not in coordsys:
            raise Exception('Invalid coordinate system: {}'.format(coordsys))
        regstrs = list(f)
    return regstrs


################
# Region parsing
################

def parse_regstr(regstr):
    """Break a region specification into pieces.
    Returns three items: rtype (str), rnums (list of floats), rprops (str)"""
    rtype, rnums, rprops = re.split('[\(\)]', regstr)  # This assumes there are no other parentheses floating around...
    rnums = [float(x) for x in re.split('\s*,\s*', rnums)]   # Could improve by only matching first set of parentheses
    return rtype, rnums, rprops


if __name__ == '__main__':
    main()

