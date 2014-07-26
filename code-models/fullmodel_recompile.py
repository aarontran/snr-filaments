"""
Quick and dirty Python script to recompile FullEfflength_mod.py
into fullmodel.so (for Python import).

Remember to run fullmodel.readfglists() before doing any computations

Aaron Tran
"""

import numpy as np
from numpy import f2py

def main():
    with open('FullEfflength_mod.f', 'r') as f:
        fsource = ''.join(list(f))
    np.f2py.compile(fsource, modulename='fullmodel', verbose=1)

if __name__ == '__main__':
    main()
