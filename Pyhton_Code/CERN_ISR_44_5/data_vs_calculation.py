#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 03:16:23 2022

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt



#From Calculations
p_N = [0.007442598,0.026021959,0.050717867,0.070385501,0.077964025,0.073472046,0.061292103,0.046477902,0.032638709,0.021516562,0.013453072,0.008041252,0.004623810,0.002570615,0.001387457,0.000729496,0.000374699,0.000188469,0.000093023]

p_N = np.array(p_N)
print(sum(p_N))
p_N =p_N/sum(p_N)

#From Data
s_N = [0.29,1.36,3.40,4.26,4.60,4.28,3.39,2.77,1.73,1.24,0.78,0.52,0.33,0.16,0.08,0.03,0.02,0.02,0.01]
err_sN = [0.05,0.14,0.36,0.38,0.34,0.35,0.27,0.21,0.14,0.09,0.06,0.06,0.07,0.05,0.03,0.02,0.02,0.02,0.01]
s_N = np.array(s_N)
err_sN = np.array(err_sN)

prob_N = s_N/sum(s_N)
err_probN = err_sN/sum(s_N)

N = np.arange(2,40,2)

print(sum((1/sum(p_N))*p_N))
print(sum(prob_N))

print(np.divide(prob_N,p_N).mean())

#writing in file
with open('/home/soumya/Softwares/pythia/pythia8307/examples/Tsallis_Calculation.txt','w') as file:
    for i in range(len(p_N)):
        file.write(str(N[i]))
        file.write('    ')
        file.write(str(p_N[i]))
        file.write('\n')
        
with open('Tsallis_probability_data.txt','w') as file:
    for i in range(len(p_N)):
        file.write(str(N[i]))
        file.write('    ')
        file.write(str(prob_N[i]))
        file.write('    ')
        file.write(str(err_probN[i]))
        file.write('\n')

plt.errorbar(N,prob_N,yerr=err_probN,fmt='o',label='CERN ISR Data')
plt.scatter(N,(1/sum(p_N))*p_N,label='Calculation from theory',color='r')
plt.xlabel('Charge Multiplicity',fontsize=20)
plt.ylabel('Probability',fontsize=20)
plt.legend()

plt.show()
