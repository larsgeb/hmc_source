import numpy as np
import random as random
import matplotlib.pyplot as plt

#============================================================
#- Read data.
#============================================================

fid=open('DATA/data.txt')
data=fid.read().strip().split('\n')
fid.close()

nt=int(data[0].split()[0])
dt=float(data[0].split()[1])
nrec=int(data[0].split()[2])

ux=np.zeros([nrec,nt])
uy=np.zeros([nrec,nt])
uz=np.zeros([nrec,nt])

cx=np.zeros([nrec,nt])
cy=np.zeros([nrec,nt])
cz=np.zeros([nrec,nt])

for i in range(nrec):

	ux[i,:]=np.array(data[1+6*i].split()).astype(np.float)
	cx[i,:]=np.array(data[2+6*i].split()).astype(np.float)
	uy[i,:]=np.array(data[3+6*i].split()).astype(np.float)
	cy[i,:]=np.array(data[4+6*i].split()).astype(np.float)
	uz[i,:]=np.array(data[5+6*i].split()).astype(np.float)
	cz[i,:]=np.array(data[6+6*i].split()).astype(np.float)

#============================================================
#- Plot data.
#============================================================

t=np.arange(0.0,nt*dt,dt)
scale=np.max([ux,uy,uz])

f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)

for i in range(nrec):
	ax1.plot(t,ux[i,:]+scale*(i+1),'k')
	ax1.plot(t,1.0/np.sqrt(cx[i,:])+scale*(i+1),'k--')
	ax1.set_title('x-component')
	ax1.set_xlabel('time [s]')

	ax2.plot(t,uy[i,:]+scale*(i+1),'k')
	ax2.plot(t,1.0/np.sqrt(cy[i,:])+scale*(i+1),'k--')
	ax2.set_title('y-component')
	ax2.set_xlabel('time [s]')

	ax3.plot(t,uz[i,:]+scale*(i+1),'k')
	ax3.plot(t,1.0/np.sqrt(cz[i,:])+scale*(i+1),'k--')
	ax3.set_title('z-component')
	ax3.set_xlabel('time [s]')

plt.show()

fid.close()

