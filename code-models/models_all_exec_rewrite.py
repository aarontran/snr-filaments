"""
Prototype for a massive redesign....
code is giving me too many headaches and is hard to work with and understand

Aaron Tran
August 2014
"""

# What if I build a Fitter class that performs table fits, etc.?
# Easiest place to globally set verbose flag. Motivation -- ugly error finding
# code...  Essentially, design with limited "global" variables...

# Principle: all these fits + errors revolve around:
# SNR, kevs, data, eps, inds for some region
# FIXED mu, eta2dict

# the rest is just iteration -- change inds, change mu/eta2dict, change SNR

# build a new Fitter() object for each dataset to be fitted.
# what are the advantages?
# * hold data all in one place (avoid long parameter lists, where most of the
#   time parameters are NOT being varied)
# * construct full/simp fitter methods for given dataset to fit.
#   this allows simpler methods for error calculation, etc?

# Fitter() object should NOT handle data presentation/formatting, ideally


#def full_fitter(snr, kevs, data, eps, mu):
#    """Wrapper function to package SNR and data dependence"""
#    return lambda **kwargs: full_fit(snr, kevs, data, eps, mu, **kwargs)

class Fitter(object):

    def __init__(self, snr, kevs, data, eps, tab):
        self.snr
        self.kevs
        self.data
        self.eps
        self.eta2dict

        # Functions
        self.vprint
        self.fitter_full
        self.fitter_simp

    def get_full_fit()                  # Return: data for table
    def get_full_err(best fit params)   # Return: data for table

    def get_simp_fit()                  # Return: data for table
    def get_simp_err(res)               # Return: data for table

    def _get_ci_errors(...)
    def _get_ci_bounds(...)
    def _get_bounds(...)
    def _one_dir_root(...)

# Code structure:
# build an SNR object
# Iterate over data:
#   modify SNR object as appropriate
#   mask data: kevs[inds], data[inds], eps[inds], tab = msk_tab(inds, tab)
#   create fitter function(snr, kevs, data, eps)
#   eta2, B0, chisqr_min = get_table_fit(fitter, tab)  # replace pars w/res
#   ci_dict = get_table_err(eta2, B0, chisqr_min, fitter, tab, **kwargs)
#   display_tab(ans, ansE)
#   latex_tab(ans, ansE)


