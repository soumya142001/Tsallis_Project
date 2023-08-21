#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 02:04:37 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma


data = np.genfromtxt('CERN_ISR_sigma_52_6.txt')

sigma = data[:,1]
d_sigma = data[:,2]

n = data[:,0]
P_n = sigma/sum(sigma)
err_Pn = d_sigma/sum(sigma)

with open('CERN_ISR_Pn_n_52_6.txt','w') as file:
    for i in range(len(P_n)):
        file.write(str(n[i]))
        file.write('    ')
        file.write(str(P_n[i]))
        file.write('    ')
        file.write(str(err_Pn[i]))
        file.write('\n')
        
file.close()