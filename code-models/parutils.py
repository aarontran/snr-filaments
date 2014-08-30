"""
Utilities for parallel computation in IPython notebooks
(don't import into a vanilla Python session, this is meant for more interactive
usage)

This is obviously not designed for high power computing, it's intended to take
advantage of all 4 cores on personal computers, nowadays.

Aaron Tran
August 2014
"""

import sys
import time

from IPython.display import clear_output
from IPython import parallel

# ============================
# Setup / boilerplate refactor
# ============================

def get_lview():
    """Obtain LoadBalancedView object for computation"""
    clients = parallel.Client()
    clients.direct_view().use_dill()  # Testing, stackoverflow.com/a/24316222
    lview = clients.load_balanced_view()
    lview.block = False
    return lview

def get_dview():
    """Obtain DirectView object for computation"""
    clients = parallel.Client()
    dview = clients[:]
    dview.block = False
    dview.use_dill()  # to serialize messy things
    return dview


# ======================================
# Stdout display for AsyncResult objects
# when running code in non-blocking mode
# ======================================

# See: http://stackoverflow.com/q/15289168
# and: http://stackoverflow.com/q/18086299
# both answered by Benjamin (Min) Ragan-Kelley

def stdout_when_done(ar, dt=1):
    """(blocking) print stdout all at once when done"""
    while not ar.ready():
        time.sleep(dt)
    ar.display_outputs()

def print_stdout(ar, dt=1, truncate=1000):
    """(blocking) print stdout while waiting for ar to finish
    Clears prior input
    Modified from http://stackoverflow.com/a/18095126
    """
    while not ar.ready():
        stdouts = ar.stdout
        if not any(stdouts):
            continue
        clear_output()
        print '-' * 30
        print "%.3fs elapsed" % ar.elapsed
        print ""
        for n, stdout in zip(xrange(len(ar.stdout)), ar.stdout):
            if stdout:
                print "[ stdout %2i ]\n%s" % (n, stdout[-truncate:])
        sys.stdout.flush()
        time.sleep(dt)


if __name__ == '__main__':
    pass
