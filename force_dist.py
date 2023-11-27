from header import *
from func import *
#########################################################
kasigma=np.array([KA,sigma_0])
t_delta,t_FF=forcedistcurve(kasigma)
fig,ax=plt.subplots()
ax.plot(t_delta,t_FF,'o-')
ax.plot(expt_data[:,0]-expt_data[0,0],expt_data[:,1],'o-')
np.savetxt(fname,np.vstack([t_delta,t_FF]))
plt.show()