# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:50:04 2014

@author: Sean
"""
import numpy
import lmfit
import matplotlib.pyplot as plt

data1 = numpy.array([35.5, 31.94, 25.34])
eps_data1 = numpy.array([1.73, .97, 1.71])
data2 = numpy.array([23.02, 17.46, 15.3])
eps_data2 = numpy.array([.35,.139, .559])
data3 = numpy.array([49.14, 42.76,29.34])
eps_data3 = numpy.array([1.5, .718, .767])
data4 = numpy.array([29, 23.9, 16.6])
eps_data4 = numpy.array([.9, .39, .45])
data5 = numpy.array([33.75, 27.2, 24.75 ])
eps_data5 = numpy.array([2.37,.62,.61])


data = data5
eps_data = eps_data5
def Width(params,data,eps_data):
    
    B0 = params['B0'].value
    eta2 = params['eta2'].value
    mu =2
    
    vd = 5e8/4.0
    b = 1.57e-3
    Cd = 2.083e19
    cm = 1.82e18
    rs = 2.96e19
    rsarc = 900.
    nu2 = 4.8359e17
    nu = numpy.array([.7*4.8359e17/2, 4.8359e17/2, 4.8359e17])
    
    eta = eta2*numpy.sqrt(nu2/cm/B0)**(1-mu)
    E = numpy.sqrt(nu/cm/B0)
    tsynch = 1.0/b/B0**2/E

    D = eta*Cd*E**mu/B0
    a = 2*D/vd/(numpy.sqrt(1+(4*D/vd**2/tsynch))-1)*4.6
    
    return (a/rs*rsarc-data)/eps_data
    
    


params = lmfit.Parameters()
params.add('B0',value = 1e-4, min = 1e-6, max = 1e-3)
params.add('eta2',value = 1, min = 0, max = 1e4)

result =lmfit.minimize(Width,params,method='leastsq',args = (data,eps_data) )

print result.chisqr

print params


final = data + result.residual

import pylab
nu = numpy.array([.7*4.8359e17/2, 4.8359e17/2, 4.8359e17])
plt.errorbar(nu,data, yerr= eps_data)

vd = 5e8/4.0
b = 1.57e-3
Cd = 2.083e19
cm = 1.82e18
rs = 2.96e19
mu = 2
nugraph = numpy.arange(.7*4.8359e17/2,4.8359e17, (1.0-.7/2.0)*4.8359e17/1000.0 )
Egraph = numpy.sqrt(nugraph/cm/B0)
B0 = params['B0'].value
eta2 = params['eta2'].value
nu2 = 4.8359e17
eta = eta2*numpy.sqrt(nu2/cm/B0)**(1-mu)
tsynch = 1.0/b/B0**2/Egraph
D = eta*Cd*Egraph**mu/B0
a = 2*D/vd/(numpy.sqrt(1+(4*D/vd**2/tsynch))-1)*4.6
    

pylab.plot(nugraph,a/rs*rsarc)

pylab.plot()