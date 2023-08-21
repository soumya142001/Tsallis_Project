#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 22:42:31 2023

@author: soumya
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import kn

N_avg = 12.08
k = 8.6
a=5.07614*10**(-3)
v0=0.368*a**3
V=40.1*a**3
g=1
g_pi=3
g_omega=1
g_eta=1
g_rho = 3
m_pi=139.5
m_eta=548.8
m_rho=775
m_omega=782
m=139.5

q1=(1/k)+(2*v0)/V
q2=(2*v0)/V

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

#def n0(u):
#    return phi(u)*V

def n0(u):
    return (phi_pi(u)+phi_eta(u)+phi_rho(u)+phi_omega(u))*V


def xi_eq(v):
    s = (q1*n0(v))/((q1-q2)*(n0(v))**2+n0(v)-N_avg)
    return s

def xi_diff(w):
    return (3+w*m*(kn(1,w*m)/kn(2,w*m)))

def xi_diff_multi(w):
    nf = g_pi*(m_pi**2)*kn(2,w*m_pi)+g_eta*(m_eta**2)*kn(2,w*m_eta)+g_rho*(m_rho**2)*kn(2,w*m_rho)+g_omega*(m_omega**2)*kn(2,w*m_omega)
    nc = 1/nf
    #print(nc)
    t_pi = g_pi*(m_pi**2)*(kn(2,w*m_pi)+(w*m_pi)*kn(1,w*m_pi)+2*kn(2,w*m_pi))
    t_rho = g_rho*(m_rho**2)*(kn(2,w*m_rho)+(w*m_rho)*kn(1,w*m_rho)+2*kn(2,w*m_rho))
    t_eta = g_eta*(m_eta**2)*(kn(2,w*m_eta)+(w*m_eta)*kn(1,w*m_eta)+2*kn(2,w*m_eta))
    t_omega = g_omega*(m_omega**2)*(kn(2,w*m_omega)+(w*m_omega)*kn(1,w*m_omega)+2*kn(2,w*m_omega))
    
    return nc*(t_pi+t_rho+t_eta+t_omega)
    
    
def q_1_diff(z):
    return q1/((xi_diff(z))**2)

def q_1_diff_multi(z):
    return q1/((xi_diff_multi(z))**2)

def q_1_eq(z):
    return q1/((xi_eq(z))**2)

def Tsallis_vdwaal(vars):
    beta_sol = vars
    #beta,q = symbols('beta, q')
    #n0 = ((g*m**3)/(2*(np.pi**2)*m*x))*kn(2,m*x)
    #chi = 3+m*x*(kn(1,m*x)/kn(2,m*x))
    
    eq1 = xi_eq(beta_sol)
    eq2 = xi_diff(beta_sol)
    
    return [eq1,eq2]


beta_inv=np.linspace(100,250,10000)
beta = np.divide(np.ones(len(beta_inv)),beta_inv)

#beta_sol =  fsolve(Tsallis_vdwaal)

#print(beta_sol)

#print(beta)

plt.plot(beta,xi_eq(beta),label=r'$\xi$'+r'($\beta$)'+' eq(33)',linewidth=4)
plt.plot(beta,xi_diff_multi(beta),label=r'$\xi$'+r'($\beta$)'+' eq(41)',linewidth=4)
#plt.plot(beta,xi_diff(beta),label="$xi$ from differentiation")
#plt.plot(beta,q_1_diff(beta),label='q-1 as a function of beta from differentiation')
#plt.plot(beta,q_1_diff_multi(beta),label='q-1 vs \u03B2 obtained by using \u03BE from eq(13)')
#plt.plot(beta,q_1_eq(beta),label='q-1 vs \u03B2 obtained by using \u03BE from eq(15)')
#plt.plot(beta,n0(beta),label='n0')
plt.xlabel(r'$\beta$',fontsize=40)
plt.ylabel(r'$\xi$',fontsize=40)
plt.ylim([-1,10])
plt.xlim([0.0061,0.00622])
plt.xticks(np.arange(6.10e-3,6.23e-3,2e-5),fontweight='bold',fontsize=20)
plt.yticks(np.arange(-1,10,1),fontweight='bold',fontsize=20)
plt.legend(frameon=False,fontsize=30)
plt.show()

print(N_avg/n0(0.006))
print(q_1_eq(0.00625))
print(q_1_diff_multi(0.00625))
print(q1)
