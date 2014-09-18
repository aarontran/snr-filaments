"""
August 5 2014

Short script to merge fragmented Tycho tables
meant for one time use (operated on pieces now stored in separate folder)
"""

import cPickle as pickle
import glob

fnames = glob.glob('Tycho_gen_2014-07-30_grid_6-100-20_vs-*.pkl')

f_orig = [f for f in fnames if 'part2' not in f]
f_pt2_1 = [f for f in fnames if 'mu-1.50' in f]
f_pt2_2 = [f for f in fnames if 'mu-2.00' in f]

f_orig.sort()
f_pt2_1.sort()
f_pt2_2.sort()

def pkldict(f):
    with open(f, 'r') as fpkl:
        d = pickle.load(fpkl)
    return d

# pt2_1 refers to mu=1.50, with subset of eta2 values
# pt2_2 refers to mu=2.00, with all eta2 values

for orig, pt2_1, pt2_2 in zip(f_orig, f_pt2_1, f_pt2_2):

    main_dict = pkldict(orig)
    main_dict[1.5].update( pkldict(pt2_1)[1.5] )
    main_dict.update( pkldict(pt2_2) )
    
    with open(orig + 'all', 'w') as fpkl:
        pickle.dump(main_dict, fpkl)
