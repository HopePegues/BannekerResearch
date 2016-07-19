# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:04:10 2016

@author: hope1707
"""

import numpy as np

#Using ratio of planet and star radius to find transit signal
#When implimented into real code the radii will be real numbers from the actuall data
plt_radius = np.array([np.random.random(100)])
star_radius = np.array([np.linspace(0.06,2,100, endpoint = True)])

percent_ratio =  plt_radius/star_radius * 100
print percent_ratio

