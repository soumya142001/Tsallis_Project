import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma

data = np.genfromtxt('/home/soumya/Tsallis_Data/Pyhton_Code/ep_data/ep-H1-data/ep-H1-data/ep_py_data/t14-epH1_py.txt')

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

print('#### Starting ep-t14 #####')
print('R_squared : '+str(R_sq))
print('N-parameter: '+str(fit_A*fit_k))
print('k-parameter:'+str(fit_k))
print('Normalisation Constant:'+str(fit_nc))

norm = (1/sum(NB(N,fit_A,fit_k,fit_nc)))

p_NB = NB(N,a,k,fit_nc)

print(norm)

x=np.linspace(1,43,1000)

#plt.plot(N,NB_fit(N,fit_A,fit_k,fit_C),label='Fitting Curve')
plt.plot(x,NB(x,fit_A,fit_k,fit_nc),label='Fitting Curve')
#plt.scatter(N,p_NB,label='Negative Binomial')
plt.errorbar(N,prob_N,yerr = err_prob_N, fmt='o',label='ep-t14 data')

plt.xlabel('Charged Multiplicity, $N_c$',fontsize=20)
plt.ylabel('Probabilty, $P_n$',fontsize=20)    
plt.legend()    
#plt.show()

figure = plt.gcf()  # get current figure
figure.set_size_inches(18,9) # set figure's size manually to your full screen (32x18)

plt.savefig('/home/soumya/Tsallis_Data/Plots/ep_plots/t14_plots/ep_t14_NB_fit.png')
