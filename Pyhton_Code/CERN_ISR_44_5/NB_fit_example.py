import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma

#Enter the .txt file or.csv file that you want to fit 
data = np.genfromtxt('Your data file name')

#Data has columns as N, error in N, Probability of N, error in Prob N

N = data[:,0]
prob_N = data[:,1]
err_N = data[:,2]
err_prob_N = data[:,3]

#Calculating Mean and std from Data 

N_avg = np.dot(N,prob_N)

N2_avg = np.dot(N**2,prob_N)

D = np.sqrt(N2_avg-(N_avg)**2)

k = (N_avg**2)/(D**2-N_avg)
a = N_avg/k

#Printing the info
print(k)
print(N_avg)

#initial guess for fit
guess=[2,7.6,2]

#The negative Binomial Function definition
def NB(x,A,K,nc):
    
     z1 = gamma(x+K)
     z2 = gamma(x+1)*gamma(K)
     z3 = A**(x)
     z4 = (1+A)**(x+K)
     
     z=nc*(z1*z3)/(z2*z4)
     
     return z
 
 
parameters,covariance = curve_fit(NB,N,prob_N,p0=guess)

#The Fitting parameters redefining 
fit_A = parameters[0]
fit_k = parameters[1]
fit_nc = parameters[2]


predict = NB(N,fit_A,fit_k,fit_nc)
actual = prob_N

#Calculating the R^2 value
corr_matrix = np.corrcoef(actual,predict)
corr = corr_matrix[0,1]
R_sq = corr**2

print('#### Starting ep-t10 #####')
print('R_squared : '+str(R_sq))
print('N-parameter: '+str(fit_A*fit_k))
print('k-parameter:'+str(fit_k))
print('Normalisation Constant:'+str(fit_nc))

norm = (1/sum(NB(N,fit_A,fit_k,fit_nc)))

p_NB = NB(N,a,k,fit_nc)

print(norm)

x=np.linspace(1,43,1000)

plt.plot(x,NB(x,fit_A,fit_k,fit_nc),label='Fitting Curve')
plt.errorbar(N,prob_N,yerr = err_prob_N, fmt='o',label='ep-t10 data')

plt.xlabel('Charged Multiplicity, $N_c$',fontsize=20)
plt.ylabel('Probabilty, $P_n$',fontsize=20)    
plt.legend()    
plt.show()


