#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:17:01 2023

@author: soumya
"""
print('## Starting PiEm525 ## ')

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import kn
from sympy import *

N_avg = 15.8469
k = 3.9667
a=5.07614*10**(-3)
v0=0.368*a**3
V=40.1*a**3
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

def n0(u):
    return (phi_pi(u)+phi_eta(u)+phi_rho(u)+phi_omega(u))*V

def xi_pi(w):
    return (3+w*m*(kn(1,w*m)/kn(2,w*m)))

def xi_eq(v):
    s = (q1*n0(v))/((q1-q2)*(n0(v))**2+n0(v)-N_avg)
    return s

def xi(w):
    nf = g_pi*(m_pi**2)*kn(2,w*m_pi)+g_eta*(m_eta**2)*kn(2,w*m_eta)+g_rho*(m_rho**2)*kn(2,w*m_rho)+g_omega*(m_omega**2)*kn(2,w*m_omega)
    nc = 1/nf
    #print(nc)
    t_pi = g_pi*(m_pi**2)*(kn(2,w*m_pi)+(w*m_pi)*kn(1,w*m_pi)+2*kn(2,w*m_pi))
    t_rho = g_rho*(m_rho**2)*(kn(2,w*m_rho)+(w*m_rho)*kn(1,w*m_rho)+2*kn(2,w*m_rho))
    t_eta = g_eta*(m_eta**2)*(kn(2,w*m_eta)+(w*m_eta)*kn(1,w*m_eta)+2*kn(2,w*m_eta))
    t_omega = g_omega*(m_omega**2)*(kn(2,w*m_omega)+(w*m_omega)*kn(1,w*m_omega)+2*kn(2,w*m_omega))
    
    return nc*(t_pi+t_rho+t_eta+t_omega)



def Tsallis_vdwaal(vars):
    beta,q = vars
    #beta,q = symbols('beta, q')
    #n0 = ((g*m**3)/(2*(np.pi**2)*m*x))*kn(2,m*x)
    #chi = 3+m*x*(kn(1,m*x)/kn(2,m*x))
    
    eq1 = n0(beta)+q*n0(beta)*xi(beta)*(n0(beta)*xi(beta)-1)-((2*v0)/V)*(n0(beta))**2-N_avg
    eq2 = q*(xi(beta))**2-2*(v0/V)-1/k
    
    return [eq1,eq2]

beta,q =  fsolve(Tsallis_vdwaal,(0.001,0.02))
#sols = nsolve(Tsallis_vdwaal,[beta,q],[1,1])

print(beta,q,beta**(-1))
print('beta : '+str(beta))
print('q-value : '+str(q))
print('temp. or beta inv. '+str(beta**(-1)))
print(Tsallis_vdwaal((beta,q)))
