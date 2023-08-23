import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sqrts = np.array([30.4,44.5,52.6,62.2])
Q = np.zeros(len(sqrts))

Q[0] = 1.0089
Q[1] = 1.0095
Q[2] = 1.0104
Q[3] = 1.0105

a,b = np.polyfit(sqrts,Q,1);


plt.scatter(sqrts,Q,s=500,marker='o') 
plt.plot(sqrts,a*sqrts+b,linewidth=4,color='r')  

plt.xlabel(r'$\sqrt{s}$ (GeV)',fontsize=40,labelpad=10)
plt.ylabel('q',fontsize=40,labelpad=20)

plt.legend(frameon=False,fontsize=30)

plt.xticks(np.arange(25,66,5),fontweight='bold',fontsize=20)
plt.yticks(np.arange(1,1.021,0.003),fontweight='bold',fontsize=20)

plt.xlim(25,65)
plt.ylim(1,1.02)

plt.show()
