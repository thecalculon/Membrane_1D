# Following are the assumptions
# R1<RI<R0
# R0>RI>R1
# The only issue is: r = 0 seems singular
# Coordinate system::
# ^ z
# |
# ------> r
from header import *
from func import *
#+begin_src python :session mempy :results output
#+end_src
#+begin_src python :session mempy :results output
#################################################################
# Lets start playing around
# R0>RI>R1
# t_R0=np.arange(1.058, 1.060, 0.0001)
# t_R0=[1.023]
# # t_R1=np.arange(1e-4, 10e-2, 1e-3)
# t_R1=[70e-4]
# # t_R0=[1.03]
# # t_R1=[8e-2]
# # FF=np.linspace(1,30,30)*1e-2
# FF=np.linspace(0.9,10,40)*1e-4
# # FF=np.array([1e-2])
# delta=np.zeros(FF.shape[0])
# t_ft=np.zeros(FF.shape[0])
# R0=np.zeros(FF.shape[0])
# RI=np.zeros(FF.shape[0])
# R1=np.zeros(FF.shape[0])
# for i,F_inp in enumerate(FF):
#    R0[i], RI[i], R1[i], delta[i], t_ft[i] = delta4force(F_inp,t_R0,t_R1,kasigma)
#    t_R0=np.arange(R0[i],R0[i]+0.04,0.001)
#    t_R1=np.arange(R1[i], R1[i]+0.03, 1e-3)
# dd=np.vstack([delta,t_ft,R0,RI,R1])
# np.savetxt("data/sigma0o001_k0o5.txt",dd)
#+end_src
#+RESULTS:
#################################################################
kasigma=np.array([KA,sigma_0])
# expt_data=np.loadtxt("data/sigma0o00076_k0o026.txt").T
# expt_data[:,0]=expt_data[:,0]-expt_data[0,0]
# cost,t_delta,t_FF=cost_func(kasigma,expt_data,Np=100)
# fig, ax = plt.subplots()
# ax.plot(expt_data[:,0],expt_data[:,1],'-s', markerfacecolor="none", label="Theory")
# ax.plot(t_delta,t_FF,'-s', markerfacecolor="none", label="Theory")
# plt.show()
# print(cost)
# t_delta,t_FF=forcedistcurve(kasigma=kasigma,Np=100)
# exit(1)
#--------------------------------------------------------------#
expt_data=np.loadtxt("data/sigma0o0001_k0o5.txt").T
expt_data=expt_data[0:20,:]
options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
initial_guess_1 = kasigma
dimensions = initial_guess_1.size
initial_guess = initial_guess_1.reshape((1, dimensions))
n_particles=2
init_pos = np.tile(initial_guess, (n_particles, 1))
optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,
		dimensions=2,options=options)
stats = optimizer.optimize(cost_func, iters=100, expt_data=expt_data,verbose=True)
print(stats)
# # #--------------------------------------------------------------#
# fname="data/sigma"+str(sigma_0).replace(".","o")+\
# 	"_k"+str(KA).replace(".","o")+".txt"
# print(fname)
# np.savetxt(fname,np.vstack([t_delta,t_FF]))
# # #--------------------------------------------------------------#
# fig, ax = plt.subplots()
# ax.plot(t_delta,t_FF,'-s', markerfacecolor="none", label="Theory")
#--------------------------------------------------------------#
# expt_data = np.loadtxt("./Ref_data/wt_expt_our.dat")
# expt_data[:,0]=expt_data[:,0]
# expt_data[:,1]=expt_data[:,1]
# radius_ev=51e-9
# start=57
# end=59
# # print(expt_data[start:end,0]/radius_ev,expt_data[start:end,1]/(sigma_0*radius_ev))
# xx=(expt_data[start:end,0]/radius_ev)
# yy=(expt_data[start:end,1]/(sigma_0*radius_ev))
# zz=np.polyfit(np.log10(xx),np.log10(yy-0.0206),1)
# print(zz[0])
# ax.plot(xx,yy,'-o',label="Expt")
#--------------------------------------------------------------#
# print((yy[1]-yy[0])/(xx[1]-xx[0]))
# zz=np.polyfit(xx, yy, 1)
# print(zz)
# p=np.poly1d(zz)
# yy=p(xx)
# ax.plot(xx,yy,'-k')
# print(zz)
# print(xx,yy)
# xx=np.array([0.007,0.20])
# yy=np.poly1d(zz)(xx)
# yy= 0.30*(xx**0.1) - 0.16
# print(xx,yy)
#--------------------------------------------------------------#
# xx=np.array([0.3,0.65])
# yy=np.array([0.2,0.6])-0.05
# ax.plot(xx,yy,'-k')
#--------------------------------------------------------------#
# xx=np.array([-2.1,-1.4])
# yy=0.25*xx-1.2
# yy=0.15*xx+expt_data[start,1]/(sigma_0*radius_ev)-0.15*xx[0]
# ax.plot(xx,yy,'-k')
#--------------------------------------------------------------#
# xx=np.array([expt_data[start+10,0]/radius_ev,expt_data[start+40,0]/radius_ev])
# yy=1.10*xx**2+expt_data[start+10,1]/(sigma_0*radius_ev)-1.10*xx[0]**2
# ax.loglog(xx,yy,'-k')
#--------------------------------------------------------------#
# # ax.set(xlim=[-0.1,2.0], ylim=[-0.01,1.0], xlabel=r"$\delta/R_v$", ylabel="F/nN")
# ax.legend()
# # plt.savefig("results/sigma1em4_k0o05.pdf")
# plt.show()
#---------------------------------------------------------------#