"""
Utilities for working w/ Python json

Aaron Tran
September 18
"""

import json
import numpy as np

import lmfit

def dump(jobj, fname, cls=None, indent=4):
    """Using my preferred defaults"""
    with open(fname, 'w') as fjson:
        json.dump(jobj, fjson, cls=cls, indent=indent)

def load(fname):
    with open(fname, 'r') as fjson:
        a = json.load(fjson)
    return a

class NumpyJSONEncoder(json.JSONEncoder):
    """Convert nd.ndarray to list for JSON encoding
    From http://stackoverflow.com/a/10367428
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class LmfitJSONEncoder(NumpyJSONEncoder):
    """Handle lmfit.Parameter() objects by saving useful properties to dict
    The selection of parameters follows lmfit.Parameter.__getstate__
    """
    def default(self, obj):
        if isinstance(obj, lmfit.Parameter):
            dobj = {}
            dobj['name'] = obj.name
            dobj['value'] = obj.value
            dobj['vary'] = obj.vary
            dobj['expr'] = obj.expr
            dobj['min'] = obj.min
            dobj['max'] = obj.max
            dobj['stderr'] = obj.stderr
            dobj['correl'] = obj.correl
            dobj['init_value'] = obj.init_value
            return dobj
        return NumpyJSONEncoder.default(self, obj)


if __name__ == '__main__':
    pass
