import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/home/vipin/.matplotlibrc')
import sys
########################################################
radius_ev=1e-5
sigma_0=1e-4
fig,ax=plt.subplots()
########## Comparison with the book #####################
# dd=np.loadtxt("data/sigma0o0001_k0o05.txt")
# dd[0,:]=dd[0,:]-dd[0,0]
# zz=np.polyfit(dd[0,:10],dd[1,:10],1)
# print(np.log10(zz[0]))
# ax.plot((dd[0,:]),(dd[1,:]-zz[1]),'s-',label=r'$K_A=0.05$')
# xx=np.array([(dd[0,0]),(dd[0,10])])
# yy=xx-0.189
# ax.plot(xx,yy,'k--')
# dd=np.loadtxt("GVB/fig12p5a_sigma0o0001_k0o05.csv")
# dd[:,0]=dd[:,0]*1e-6
# dd[:,1]=dd[:,1]*1e-9
# ax.plot((-dd[:,0]/radius_ev),(dd[:,1]/(radius_ev*sigma_0)),'.-k')
# # #------------------------------------------------------------#
# dd=np.loadtxt("data/sigma0o0001_k0o25.txt")
# dd[0,:]=dd[0,:]-dd[0,0]
# dd[1,:]=dd[1,:]-dd[1,0]
# zz=np.polyfit(dd[0,:10],dd[1,:10],1)
# print(np.log10(zz[0]))
# ax.plot(dd[0,:],dd[1,:],'s-',label=r'$K_A=0.25$')
# dd=np.loadtxt("GVB/fig12p5a_sigma0o0001_k0o25.csv",delimiter=',')
# dd[:,0]=dd[:,0]*1e-6
# dd[:,1]=dd[:,1]*1e-9
# ax.plot(-dd[:,0]/radius_ev,dd[:,1]/(radius_ev*sigma_0),'.-k')
# # #------------------------------------------------------------#
# dd=np.loadtxt("data/sigma0o0001_k0o5.txt")
# dd[0,:]=dd[0,:]-6.96380422e-02
# dd[1,:]=dd[1,:]-4.95395617e-05
# zz=np.polyfit(dd[0,:10],dd[1,:10],1)
# print(np.log10(zz[0]))
# ax.plot(dd[0,:],dd[1,:],'s-',label=r'$K_A=0.5$')
# dd=np.loadtxt("GVB/fig12p5a_sigma0o0001_k0o5.csv",delimiter=',')
# dd[:,0]=dd[:,0]*1e-6
# dd[:,1]=dd[:,1]*1e-9
# print(dd[:,0])
# ax.plot(-(dd[:,0])/radius_ev,dd[:,1]/(radius_ev*sigma_0),'.-k')
#------------------------------------------------------------#
dd=np.loadtxt("GVB/GUV.csv",delimiter=',')
radius_ev=11.4e-6
xx=np.log10(-dd[0:2,0]*1e-6/radius_ev)
yy=np.log10(dd[0:2,1]*1e-9/radius_ev)
ax.plot(xx,yy,'o-')
zz=np.polyfit(xx,yy,1)
print(zz)
#------------------------------------------------------------#
# ax.set(ylim=[0,10])
ax.grid(True)
ax.legend()
plt.show()