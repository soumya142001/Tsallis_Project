#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 03:06:12 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma

data = np.genfromtxt('/home/soumya/Tsallis_Data/Pyhton_Code/had_nucl_code/had_nucl_data/had_nucl_data/piNe_30.txt')

N = data[:,0]
prob_N = data[:,1]
err_N = data[:,2]
err_prob_N = data[:,3]

#print(data)

N_avg = np.dot(N,prob_N)

N2_avg = np.dot(N**2,prob_N)

D = np.sqrt(N2_avg-(N_avg)**2)

k = (N_avg**2)/(D**2-N_avg)
a = N_avg/k

print(k)
print(N_avg)

guess=[2,7.6,2]

def Tsallis(x,beta,q,NC):
     
    w1 = (1-beta*(q-1)*x)
    w2 = NC*w1**(q/(1-q))
    
    return w2

def Tsallis_modified(x,q,NC):
     
    w1 = (1+(q-1)*x)
    w2 = NC*w1**(-q/(q-1))
    
    return w2


def NB(x,A,K,nc):
    
     z1 = gamma(x+K)
     z2 = gamma(x+1)*gamma(K)
     z3 = A**(x)
     z4 = (1+A)**(x+K)
     
     z=nc*(z1*z3)/(z2*z4)
     
     return z
 
 
parameters,covariance = curve_fit(NB,N,prob_N,p0=guess)


fit_A = parameters[0]
fit_k = parameters[1]
fit_nc = parameters[2]


predict = NB(N,fit_A,fit_k,fit_nc)
actual = prob_N

corr_matrix = np.corrcoef(actual,predict)
corr = corr_matrix[0,1]
R_sq = corr**2

print('R_squared : '+str(R_sq))
print('N-parameter: '+str(fit_A*fit_k))
print('k-parameter:'+str(fit_k))
print('Normalisation Constant:'+str(fit_nc))

norm = (1/sum(NB(N,fit_A,fit_k,fit_nc)))

p_NB = NB(N,a,k,fit_nc)

print(norm)

s='$\overline{N}$ from Data: '+str(N_avg)+'\n'+'k from Data: '+str(k)+'\n'+'$\overline{N}$ from Fit: '+str(fit_A*fit_k)+'\n'+'k from Data: '+str(fit_k)+'\n'+'$R^2$ of the fit: '+str(R_sq)

x=np.linspace(1,60,1000)

plt.text(30,0.08,s,fontsize=15,bbox=dict(facecolor='orange', alpha=0.5))
#plt.plot(N,NB_fit(N,fit_A,fit_k,fit_C),label='Fitting Curve')
plt.plot(x,NB(x,fit_A,fit_k,fit_nc),label='Fitting Curve')
#plt.scatter(N,p_NB,label='Negative Binomial')
plt.errorbar(N,prob_N,yerr = err_prob_N, fmt='o',label='piNe_30 data')

plt.xlabel('Charged Multiplicity, $N_c$',fontsize=20)
plt.ylabel('Probabilty, $P_n$',fontsize=20)    
plt.legend()    
plt.show()
