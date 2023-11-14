# Following are the assumptions
# R1<RI<R0
# R0>RI>R1
# The only issue is: r = 0 seems singular
# Coordinate system::
# ^ z
# |
# ------> r
import numpy as np
import matplotlib.pyplot as plt
from func import *
from params import *
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
#    R0[i], RI[i], R1[i], delta[i], t_ft[i] = force_dist(F_inp,t_R0,t_R1)
#    t_R0=np.arange(R0[i],R0[i]+0.04,0.001)
#    t_R1=np.arange(R1[i], R1[i]+0.03, 1e-3)
# dd=np.vstack([delta,t_ft,R0,RI,R1])
# np.savetxt("data/sigma0o001_k0o5.txt",dd)
#+end_src
#+RESULTS:
#################################################################
Np=50
t_R1=np.linspace(0.018,0.2,Np)
t_R0=np.arange(0.8,1.10,0.001)
t_delta=np.zeros(Np)
t_FF=np.zeros(Np)
for i,R1 in enumerate(t_R1):
   for R0 in t_R0:
      RI=getRI(R0,R1)
      V1, V2, V3 =  Volume(R0, RI, R1, Np=4096)
      Vol=V1+V2+V3
      cdt3 = abs(Vol - VolT) < 0.01*VolT;
      if cdt3:
         break;
   t_delta[i]=delta(R0,RI,R1)/Rv
   t_FF[i]=force_all(R0, RI, R1)[0]/(sigma_0*Rv)
   t_R0=np.arange(R0,R0+0.01,0.001)
   print(t_delta[i],t_FF[i],R0)
np.savetxt("FF.txt",t_FF)
np.savetxt("delta.txt",t_delta)
plt.plot(t_delta,t_FF,'o-')
plt.show()