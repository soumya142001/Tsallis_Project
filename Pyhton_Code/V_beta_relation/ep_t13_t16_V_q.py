#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:24:14 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import kn
from sympy import *

#N_avg = 2.4544
#k = 3.283
a=5.07614*10**(-3)
v0=0.368*a**3
#V=40.1*a**3
g=1
g_pi=3
g_eta=1
g_omega=1
g_rho=3
m_pi=139.5
m_eta=548.8
m_rho=770
m_omega=782
m=139.5

gamma_mono = 1.66
gamma_di = 1.4
#q1=(1/k)+(2*v0)/V
#q2=(2*v0)/V

V1=30
V2=150
V = np.linspace(V1,V2,1200);
Nbar=np.array([2.65,5.32,7.63,8.78])
K = np.array([3.13,3.96,6.70,12.73])

V= V*(a**3)

def phi(x):
    t = (g/(2*x*np.pi**2))*(m**2)*kn(2,m*x)
    return t

def phi_pi(x):
    t = (g_pi/(2*x*np.pi**2))*(m_pi**2)*kn(2,m_pi*x)
    return t

def phi_eta(x):
    t = (g_eta/(2*x*np.pi**2))*(m_eta**2)*kn(2,m_eta*x)
    return t

def phi_rho(x):
    t = (g_rho/(2*x*np.pi**2))*(m_rho**2)*kn(2,m_rho*x)
    return t

def phi_omega(x):
    t = (g_omega/(2*x*np.pi**2))*(m_omega**2)*kn(2,m_omega*x)
    return t

def n0(u,V):
    return (phi_pi(u)+phi_eta(u)+phi_rho(u)+phi_omega(u))*V

def xi_pi(w):
    return (3+w*m*(kn(1,w*m)/kn(2,w*m)))

#def xi_eq(v):
#    s = (q1*n0(v))/((q1-q2)*(n0(v))**2+n0(v)-N_avg)
#    return s

def xi(w):
    nf = g_pi*(m_pi**2)*kn(2,w*m_pi)+g_eta*(m_eta**2)*kn(2,w*m_eta)+g_rho*(m_rho**2)*kn(2,w*m_rho)+g_omega*(m_omega**2)*kn(2,w*m_omega)
    nc = 1/nf
    #print(nc)
    t_pi = g_pi*(m_pi**2)*(kn(2,w*m_pi)+(w*m_pi)*kn(1,w*m_pi)+2*kn(2,w*m_pi))
    t_rho = g_rho*(m_rho**2)*(kn(2,w*m_rho)+(w*m_rho)*kn(1,w*m_rho)+2*kn(2,w*m_rho))
    t_eta = g_eta*(m_eta**2)*(kn(2,w*m_eta)+(w*m_eta)*kn(1,w*m_eta)+2*kn(2,w*m_eta))
    t_omega = g_omega*(m_omega**2)*(kn(2,w*m_omega)+(w*m_omega)*kn(1,w*m_omega)+2*kn(2,w*m_omega))

    return nc*(t_pi+t_rho+t_eta+t_omega)

def Ideal_Temp(V,gamma,k):
    z = k/(V**(gamma-1))
    return z

    

temp = np.zeros((len(Nbar),len(V)))
Q = np.zeros((len(Nbar),len(V)))

for j in range(len(Nbar)):
    N_avg = Nbar[j]
    k = K[j]
    for i in range(len(V)):
        def Tsallis_vdwaal(vars):
            beta,q = vars
            eq1 = n0(beta,V[i])+q*n0(beta,V[i])*xi(beta)*(n0(beta,V[i])*xi(beta)-1)-((2*v0)/V[i])*(n0(beta,V[i]))**2-N_avg
            eq2 = q*(xi(beta))**2-2*(v0/V[i])-1/k
        
            return [eq1,eq2]
    
        beta,q =  fsolve(Tsallis_vdwaal,(0.008,0.02))
        temp[j][i] = beta**(-1)
        Q[j][i] = q

plt.plot(V/(a**3),Q[0][:],label=r'$\langle Q^2 \rangle=13.5$ GeV$^2$',linewidth=4)
plt.plot(V/(a**3),Q[1][:],label=r'$\langle Q^2 \rangle=27.6$ GeV$^2$',linewidth=4)
plt.plot(V/(a**3),Q[2][:],label=r'$\langle Q^2 \rangle=54.9$ GeV$^2$',linewidth=4)
plt.plot(V/(a**3),Q[3][:],label=r'$\langle Q^2 \rangle=374.2$ GeV$^2$',linewidth=4)

plt.xlabel('V (fm$^3$)',fontsize=40)
plt.ylabel('q-1',fontsize=40)
plt.xticks(np.arange(20,161,20),fontweight='bold',fontsize=20)
plt.yticks(np.arange(0,0.051,0.005),fontweight='bold',fontsize=20)
plt.legend(frameon=False,fontsize=20)
plt.show()
