import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/home/vipin/.matplotlibrc')
import sys
########################################################
radius_ev=1e-9
########################################################
fig,ax=plt.subplots()
dd = np.loadtxt("./Ref_data/WT/0616ev402.txt")
ind=dd[:,0]>=0
xx=dd[ind,0]/radius_ev
yy=dd[ind,1]/radius_ev
# print(xx,yy)
ax.plot(np.log10(xx[1:]),
	np.log10(yy[1:]),'o-')
ax.plot(np.log10(xx[3:15]),np.log10(yy[3:15]),'-')

zz=np.polyfit(np.log10(xx[3:15]),np.log10(yy[3:15]),1)
print(zz)


dd = np.loadtxt("./Ref_data/CD63/1219ev906.txt")
ind=dd[:,0]>=0
xx=dd[ind,0]/radius_ev
yy=dd[ind,1]/radius_ev
# print(xx,yy)
ax.plot(np.log10(xx),
	np.log10(yy),'o-')
XX=np.array([0.1,1.5])
YY=1*XX-1.90
plt.plot(XX,YY,'--')

plt.xlabel(r"$\log (\delta\times 1e9)$")
plt.ylabel(r"$\log (F\times 1e9)$")

ax.grid(True)
plt.show()