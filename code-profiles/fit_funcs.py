
# coding: utf-8

## Profile fitting functions

# For use in scripts to, uh, fit thin synchrotron profiles!
# 
# EXPORT this notebook to .py, for import, in 'File' - 'Download as' - 'Python (.py)'.  Therefore don't use iPython magics in any cells, or delete them when done.

# In[1]:

import numpy as np


### Ressler: two-exponential + Gaussian

# First, we need a profile function.  _Ressler et al._ [2014, _in press_] used a piecewise function to fit the upstream and downstream sides of the rim profiles:
# 
# \begin{align*}
#     h_{\mathrm{up}}(x) &= A_u \exp \left(\frac{x_0-x}{w_u}\right) + C_u \\
#     h_{\mathrm{down}}(x) &= A_d \exp \left(\frac{x - x_0}{w_d}\right)
#                           + B \exp \left( \frac{-(x-x_1)^2}{2\pi \sigma^2} \right)
#                           + C_d
# \end{align*}
# 
# Some outstanding questions -- how do you decide the location of up/down-stream?  Is this part of the fitting procedure?
# The peak intensity may be located at either $x_0$ or $x_1$ -- okay, seems physically reasonable given projection effects.  But how to decide where to split the fit.

# In[2]:

def ressler(x, x0, x1, Cu, Cd, wu, wd, sigma, Au, B):
    """Fit function used by Satoru for SN 1006 profiles, see Ressler et al. [ApJ, in press]
    
    Inputs (all floats):
        x = radial position (arcsec.)
        Au, x0, wu, Cu = parameters for upstream function (x0 also relevant to downstream func)
        xs = location of upstream/downstream divide
        Ad, wd, B, x1, sigma, Cd = parameters for downstream function
    Output (float): value of function at radial position x
    """
    Ad = Au + Cu - Cd - B * np.exp( -(x0-x1)**2 / (2*np.pi*sigma**2))
    
    def ressler_help(x, x0, x1, Cu, Cd, wu, wd, sigma, Au, B, Ad):
        if x > x0:  # Upstream, ahead of shock
            return Au * np.exp((x0-x)/wu) + Cu
        else:
            return Ad * np.exp((x-x0)/wd) + B * np.exp( -(x-x1)**2 / (2*np.pi*sigma**2)) + Cd
    
    ressler_vec = np.vectorize(ressler_help)
    return ressler_vec(x, x0, x1, Cu, Cd, wu, wd, sigma, Au, B, Ad)


### Two exponential function with free split, $x_s$

# I first try a continuous piecewise two exponential function of form:
# 
# \begin{align*}
#     h_{\mathrm{up}}(x) &= A_u \exp \left(\frac{x_0-x}{w_u}\right) + C_u \\
#     h_{\mathrm{down}}(x) &= A_d \exp \left(\frac{x - x_0}{w_d}\right) + C_d
# \end{align*}
# 
# where the upstream/downstream split occurs at $x=x_s$ and $x_s$ is used as a fit parameter instead of $A_u$.  Once $x_s$ is determined, $A_u$ is constrained to be:
# 
# \begin{equation*}
#     A_u = \exp\left(\frac{x_s - x_0}{w_u}\right) \left[ A_d \exp\left(\frac{x_s - x_0}{w_d}\right) + C_d - C_u \right]
# \end{equation*}
# 
# Alternately we may constrain $A_d$ and use $A_u$ as a fit parameter:
# 
# \begin{equation*}
#     A_d = \exp\left(\frac{x_0 - x_s}{w_d}\right) \left[ A_u \exp\left(\frac{x_0 - x_s}{w_u}\right) + C_u - C_d \right]
# \end{equation*}

# In[3]:

def two_exp(x, xs, x0, Cu, Cd, wu, wd, Au):
    """Two exponential profile function, vectorized in x
    downstream amplitude Ad set by the split xs + continuity requirement
    
    Inputs (all floats):
        x = radial position (arcsec.)
        xs = location of upstream/downstream divide (and, piecewise function split)
        x0 = exponential split centering... sort of
        Cu, Cd = upstream/downstream offsets, backgrounds
        wu, wd = upstream/downstream e-folding lengthscales
        Au = upstream amplitude.  With other parameters, determines downstream amplitude Ad by continuity
    Output (float): value of function at radial position x
    """
    # Au = np.exp((xs-x0)/wu) * ( Ad * np.exp((xs-x0)/wd) + Cd - Cu )
    Ad = np.exp((x0-xs)/wd) * ( Au * np.exp((x0-xs)/wu) + Cu - Cd )
    
    def f_help(x, xs, x0, Cu, Cd, wu, wd, Au, Ad):
        if x > xs:  # Upstream, ahead of shock
            return Au * np.exp((x0-x)/wu) + Cu
        else:
            return Ad * np.exp((x-x0)/wd) + Cd
    
    f_vec = np.vectorize(f_help)
    return f_vec(x, xs, x0, Cu, Cd, wu, wd, Au, Ad)


### Two exponential function, decoupled (split $x_0$ into $x_u$ and $x_d$)

# I could also let one more parameter run free in the two exponential function:
# 
# \begin{align*}
#     h_{\mathrm{up}}(x) &= A_u \exp \left(\frac{x_u-x}{w_u}\right) + C_u \\
#     h_{\mathrm{down}}(x) &= A_d \exp \left(\frac{x - x_d}{w_d}\right) + C_d
# \end{align*}
# 
# Again, $x=x_s$ is a fit parameter for the upstream/downstream split, which can be used in place of $A_u$ or $A_d$.
# Because the upstream side is generally easier to fit / the data is better, it's probably better to pin $A_d$ and let $A_u$ be a free fit parameter.
# 
# \begin{equation*}
#     A_u = \exp\left(\frac{x_s - x_u}{w_u}\right) \left[ A_d \exp\left(\frac{x_s - x_d}{w_d}\right) + C_d - C_u \right]
# \end{equation*}
# 
# \begin{equation*}
#     A_d = \exp\left(\frac{x_d - x_s}{w_d}\right) \left[ A_u \exp\left(\frac{x_u - x_s}{w_u}\right) + C_u - C_d \right]
# \end{equation*}

# In[4]:

def two_exp_decoupled(x, xs, xu, xd, Cu, Cd, wu, wd, Au):
    """Two exponential profile function, with extra free parameter, vectorized in x
    downstream amplitude Ad set by the split xs + continuity requirement
    
    Inputs (all floats):
        x = radial position (arcsec.)
        xs = location of upstream/downstream divide (and, piecewise function split)
        xu, xd = upstream/downstream exponential locations
        Cu, Cd = upstream/downstream offsets, backgrounds
        wu, wd = upstream/downstream e-folding lengthscales
        Au = upstream amplitude.  With other parameters, determines downstream amplitude Ad by continuity
    Output (float): value of function at radial position x
    """
    # Au = np.exp((xs-xu)/wu) * ( Ad * np.exp((xs-xd)/wd) + Cd - Cu )
    Ad = np.exp((xd-xs)/wd) * ( Au * np.exp((xu-xs)/wu) + Cu - Cd )
    
    def f_help(x, xs, xu, xd, Cu, Cd, wu, wd, Au, Ad):
        if x > xs:  # Upstream, ahead of shock
            return Au * np.exp((xu-x)/wu) + Cu
        else:
            return Ad * np.exp((x-xd)/wd) + Cd
    
    f_vec = np.vectorize(f_help)
    return f_vec(x, xs, xu, xd, Cu, Cd, wu, wd, Au, Ad)


### Two exponential function, simplified; split fixed at $x_0$

# Now instead of using a free split location $x_s$, try throwing that away and using $x_0$ as the split location (following what Satoru did?):
# 
# \begin{align*}
#     h_{\mathrm{up}}(x) &= A_u \exp \left(\frac{x_0-x}{w_u}\right) + C_u \\
#     h_{\mathrm{down}}(x) &= A_d \exp \left(\frac{x - x_0}{w_d}\right) + C_d
# \end{align*}
# 
# where the upstream/downstream split occurs at $x=x_0$.  Continuity then allows us constrain one of $\{A_u, A_d, C_u, C_d\}$.
# Consistent with above, let's constrain $A_d$ and let $A_u$ run free (try changing this up as well).
# 
# \begin{align*}
#     A_d &= A_u + (C_u - C_d) \\
#     A_u &= A_d + (C_d - C_u)
# \end{align*}

# In[5]:

def two_exp_simp(x, x0, Cu, Cd, wu, wd, Au):
    """Simplified two exponential profile function, vectorized in x
    downstream amplitude Ad set by the split x0 + continuity requirement
    
    Inputs (all floats):
        x = radial position (arcsec.)
        x0 = exponential split centering, and location of piecewise function split
        Cu, Cd = upstream/downstream offsets, backgrounds
        wu, wd = upstream/downstream e-folding lengthscales
        Au = upstream amplitude.  With other parameters, determines downstream amplitude Ad by continuity
    Output (float): value of function at radial position x
    """
    # Au = Ad + Cd - Cu
    Ad = Au + Cu - Cd
    
    def f_help(x, x0, Cu, Cd, wu, wd, Au, Ad):
        if x > x0:  # Upstream, ahead of shock
            return Au * np.exp((x0-x)/wu) + Cu
        else:
            return Ad * np.exp((x-x0)/wd) + Cd
    
    f_vec = np.vectorize(f_help)
    return f_vec(x, x0, Cu, Cd, wu, wd, Au, Ad)


### Gaussian and a ramp function (just a test...)

# Another simple model, defined as:
# \begin{equation}
#     h(x) = A_0 \exp\left( \frac{-(x-x_0)^2}{2 \sigma^2} \right) + \mathrm{ramp}(x)
# \end{equation}
# where $\mathrm{ramp}(x)$ is defined as
# \begin{equation}
#     \mathrm{ramp}(x) =
#     \begin{cases}
#         C_d &\text{if }\; x < x_0 - \sigma \\
#         C_u + \frac{C_u - C_d}{2\sigma} (x - (x_0 + \sigma)) &\text{if }\; x_0 - \sigma \leq x \leq x_0 + \sigma \\
#         C_u &\text{if }\; x_0 - \sigma < x
#     \end{cases}
# \end{equation}
# 
# Of course, the ramp breakpoints may be something other than $\sigma$.

# In[6]:

def gaussramp(x, x0, Cu, Cd, sigma, a0):
    """Ramp + Gaussian function, ramp breakpoints at Gaussian 1-sigma (could free for better fit)
    Vectorized in x
    Inputs (all floats):
        x = radial position (arcsec.)
        x0 = center of Gaussian and ramp function (midway between breakpoints) (arcsec)
        Cu, Cd = upstream, downstream offsets
        sigma = Gaussian standard deviation
        a0 = Gaussian amplitude
    Output (float): value of function at radial position x
    """
    
    def ramp_help(x, x0, Cu, Cd, sigma):
        if x > x0 + abs(sigma):
            return Cu
        elif x < x0 - abs(sigma):
            return Cd
        else:
            return Cu + (x - x0 - abs(sigma)) * (Cu - Cd) / (2*abs(sigma))
    
    ramp = np.vectorize(ramp_help)
    return a0 * np.exp( -(x-x0)**2 / (2*sigma**2) ) + ramp(x, x0, Cu, Cd, sigma)

