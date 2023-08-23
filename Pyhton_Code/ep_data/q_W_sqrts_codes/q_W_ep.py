import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

eta = np.array(['1<'+r'$\eta^{*}$'+'<5','1<'+r'$\eta^{*}$'+'<4','1<'+r'4\eta^{*}$'+'<3','1<'+r'$\eta^{*}$'+'<2'])
W = np.array([96.9,132,166.8,201.9])

Q = np.zeros((len(eta),len(W)))

Q[0][0] = 1.0057
Q[1][0] = 1.0080
Q[2][0] = 1.0152
Q[3][0] = 1.0228

Q[0][1] = 1.0063
Q[1][1] = 1.0096
Q[2][1] = 1.0187
Q[3][1] = 1.0243

Q[0][2] = 1.0069
Q[1][2] = 1.0112
Q[2][2] = 1.0197
Q[3][2] = 1.0243

Q[0][3] = 1.0068
Q[1][3] = 1.0119
Q[2][3] = 1.0193
Q[3][3] = 1.0239

a0,b0 = np.polyfit(W,Q[0],1);
a1,b1 = np.polyfit(W,Q[1],1);
a2,b2 = np.polyfit(W,Q[2],1);
a3,b3 = np.polyfit(W,Q[3],1);


plt.scatter(W,Q[0],s=300,label='1<'+r'$\eta^{*}$'+'<5',marker='o') 
plt.scatter(W,Q[1],s=300,label='1<'+r'$\eta^{*}$'+'<4',marker='^') 
plt.scatter(W,Q[2],s=300,label='1<'+r'$\eta^{*}$'+'<3',marker='s') 
plt.scatter(W,Q[3],s=500,label='1<'+r'$\eta^{*}$'+'<2',marker='*') 
plt.plot(W,a0*W+b0,linewidth=4) 
plt.plot(W,a1*W+b1,linewidth=4) 
plt.plot(W,a2*W+b2,linewidth=4) 
plt.plot(W,a3*W+b3,linewidth=4) 

plt.xlabel(r'$\langle W \rangle$ (GeV)',fontsize=40,labelpad=10)
plt.ylabel('q',fontsize=40,labelpad=20)

plt.legend(frameon=False,fontsize=30,loc='right')

plt.xticks(np.arange(90,231,20),fontweight='bold',fontsize=20)
plt.yticks(np.arange(1,1.03,0.003),fontweight='bold',fontsize=20)

plt.xlim(90,250)
plt.ylim(1,1.03)

plt.show()
