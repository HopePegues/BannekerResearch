# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:12:19 2016

@author: hope1707
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate

#Using Morton&Swift's data on planet periods to make its cdf
xdata,ydata = np.loadtxt('logp_distribution_tabulated.txt', unpack=True)

cmp = np.zeros_like(xdata)

for i,ix in enumerate(xdata):
    hi = np.trapz(ydata[:i+1],xdata[:i+1])
    cmp[i] = hi
cmp = cmp/np.trapz(ydata,xdata)

spline = interpolate.splrep(cmp,xdata)
rad = interpolate.splev(np.random.uniform(0,1,10000),spline)

'''
I believe this next line does the same as everything I did above
When graphed on top of each other the two overlap 
but i am not sure if this fixes the fact that the xdata is in log scale
'''
other = integrate.cumtrapz(ydata,xdata,initial = 0)

y = np.linspace(0,1,len(xdata))
plt.plot(xdata,cmp)
plt.plot(xdata,other)
#plt.scatter(xdata,ydata)
plt.show()