import numpy as np
import random as random
import matplotlib.pyplot as plt

#============================================================
#- Setup.
#============================================================

#- Dimensions of interest.
dim_1=6
dim_2=7

#============================================================
#- Read samples and plot trajectory..
#============================================================

fid=open('OUTPUT/trajectory.txt')

dummy=fid.read().strip().split()
dimension=int(dummy[0])
#iterations=int(dummy[1])
iterations=(len(dummy)-2)/dimension

x=np.zeros(iterations)
y=np.zeros(iterations)

for i in range(iterations):

	x[i]=float(dummy[2+dim_1+i*dimension])
	y[i]=float(dummy[2+dim_2+i*dimension])

plt.plot(x,y,'k')
plt.plot(x,y,'ro')
plt.xlabel('m'+str(dim_1+1))
plt.ylabel('m'+str(dim_2+1))
plt.savefig('OUTPUT/Hamiltonian_trajectory.png')
#plt.show()

plt.close()
fid.close()