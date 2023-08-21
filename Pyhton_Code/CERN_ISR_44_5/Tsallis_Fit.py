#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 15:19:31 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma

data = np.genfromtxt('Tsallis_probability_data.txt')
data1 = np.genfromtxt('/home/soumya/Softwares/pythia/pythia8307/examples/Tsallis_Calculation.txt')

N = data[:,0]
prob_N = data[:,1]
err_prob_N = data[:,2]

k = 8.6
N_avg = 12.08
V=40.1
v0 = 0.368

print(np.divide(v0*N,(np.ones(len(N))-prob_N)))

a = N_avg/k

b = (1/k)+(2*v0)/V

guess=[1,3,2]

def Tsallis(x,beta,q,NC):
     
    w1 = (1-beta*(q-1)*x)
    w2 = NC*w1**(q/(1-q))
    
    return w2

def Tsallis_modified(x,q,NC):
     
    w1 = (1+(q-1)*x)
    w2 = NC*w1**(-q/(q-1))
    
    return w2

def Hagedron(x,beta,q,NC):
     
    T0=185.27
    y = T0*x**(11/5)
    w1 = (1-beta*(q-1)*y)
    w2 = NC*w1**(q/(1-q))
    
    return w2

def func(x,A,B,n,C):
        
     l = N_avg/(((V*n)**2)*(A*B-(2*v0)/V))
     
    
     y1 = gamma(x+l)
     y2 = gamma(x+1)*gamma(l)
     y3 = (a**(x))/((1+a)**(x+l))
     #y4 = (1+m)**(x+l)
     
     y = (C)*(y1*y3)/(y2)
    
     return y

def NB(x,A,K,nc):
    
     z1 = gamma(x+K)
     z2 = gamma(x+1)*gamma(K)
     z3 = A**(x)
     z4 = (1+A)**(x+K)
     
     z=nc*(z1*z3)/(z2*z4)
     
     return z
 
def NB_fit(x,A,K,C):
    
     z1 = gamma(x+K)
     z2 = gamma(x+1)*gamma(K)
     z3 = A**(x)
     z4 = (1+A)**(x+K)
     
     z=(C*(z1*z3))/(z2*z4)
     
     return z
 
 


    



#parameters,covariance = curve_fit(Hagedron,N,prob_N,p0=guess)
parameters,covariance = curve_fit(NB,N,prob_N,p0=guess)

# fit_A = parameters[0]
# fit_B = parameters[1]
# fit_n = parameters[2]
# fit_C = parameters[3]


# fit_beta = parameters[0]
# fit_q = parameters[1]
# fit_NC = parameters[2]

fit_A = parameters[0]
fit_k = parameters[1]
fit_nc = parameters[2]
#fit_C = parameters[2]


predict = NB(N,fit_A,fit_k,fit_nc)
actual = prob_N

corr_matrix = np.corrcoef(actual,predict)
corr = corr_matrix[0,1]
R_sq = corr**2


# print(R_sq)
print(fit_A)
print(fit_k)
print(fit_nc)
#print(fit_C)
# #print(fit_n)
#print(fit_C)

# print('R_squared : '+str(R_sq))
# print('Beta : ' +str(fit_beta))
# print('q-parameter : ' +str(fit_q))
# print('Normalisation COnstant : '+str(fit_NC))

print('R_squared : '+str(R_sq))
print('N-parameter: '+str(fit_A*fit_k))
print('k-parameter:'+str(fit_k))
#print('Normalisation Constant : '+str(fit_C))

norm = (1/sum(NB(N,fit_A,fit_k,fit_nc)))

p_NB = norm*NB(N,a,k,fit_nc)

print(norm)

#print(1+fit_A/fit_B)


with open('/home/soumya/Softwares/pythia/pythia8307/examples/Tsallis_Calculation_python.txt','w') as file:
    for i in range(len(p_NB)):
        file.write(str(N[i]))
        file.write('    ')
        file.write(str(p_NB[i]))
        file.write('\n')
        
file.close()

x=np.linspace(0,43,1000)

#plt.plot(N,NB_fit(N,fit_A,fit_k,fit_C),label='Fitting Curve')
plt.plot(x,NB(x,fit_A,fit_k,fit_nc),label='Fitting Curve')
#plt.scatter(N,p_NB,label='Negative Binomial')
plt.errorbar(N,prob_N,yerr = err_prob_N, fmt='o',label='CERN ISR Data')

plt.xlabel('Charged Multiplicity, $N_c$',fontsize=20)
plt.ylabel('Probabilty, $P_n$',fontsize=20)    
plt.legend()    
plt.show()
