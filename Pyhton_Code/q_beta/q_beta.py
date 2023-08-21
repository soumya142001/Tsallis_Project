#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 21:00:14 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

q_isr = [1.00894,1.00949,1.0104,1.0105]
sqrt_s_isr = [30.4,44.5,52.6,62.2]
beta_isr = [158.0968,161.2445,161.4110,163.1975]

m_isr,b_isr = np.polyfit(sqrt_s_isr,q_isr,1)

print([m_isr,b_isr])

s='Slope: '+str(m_isr)+'\n'+'Intercept: '+str(b_isr)

q_piem = [1.0183,1.0194,1.0133,1.0191]
sqrt_s_piem = [200,340,50,525]
beta_piem = [150.3348,151.5221,146.4367,157.6432]

q_pem = [1.0227,1.0124,1.0212,1.0194,1.0142,1.0192]
sqrt_s_pem = [200,27,300,400,67,800]

q_ep = [1.0228,1.0152,1.0080,1.0057]
Q2_ep = [13.5,27.6,55.0,385.3]

q_isr = np.array(q_isr)
sqrt_s_isr = np.array(sqrt_s_isr)
beta_isr = np.array(beta_isr)

q_piem = np.array(q_piem)
sqrt_s_piem = np.array(sqrt_s_piem)
beta_piem = np.array(beta_piem)

plt.text(35,1.01,s,fontsize=15)
# plt.scatter(beta_isr,q_isr,label='beta vs q for CERN ISR data')
plt.scatter(sqrt_s_isr,q_isr,label='sqrt{s} vs q for CERN ISR data')
plt.plot(sqrt_s_isr,m_isr*sqrt_s_isr+b_isr,label='linear fot to CERN ISR data ')
# plt.xlabel('beta')
# plt.ylabel('q-value')
# plt.title('CERN ISR data')

#plt.scatter(beta_piem,q_piem,label='beta vs q for PiEm data')
#plt.scatter(sqrt_s_piem,q_piem,label='sqrt{s} vs q for PiEm data')
#plt.scatter(sqrt_s_pem,q_pem,label='sqrt{s} vs q for PiEm data')
#plt.scatter(Q2_ep,q_ep,label='q vs $Q^2$ for ep-t1 data')
#plt.xlabel('$Q^2$')
plt.xlabel('$\sqrt{s}$')
plt.ylabel('q-value')
#plt.title('ep-t1 data')
#plt.title('PiEm data')
plt.title('CERN ISR data')

#plt.savefig('/home/soumya/Tsallis_Data/Plots/q_beta/beta_v_q_piem.png')
