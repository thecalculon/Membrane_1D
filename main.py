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
t_R0=np.arange(1, 1.05, 0.01)
t_R1=np.arange(1e-4, 0.02, 1e-3)
# FF=np.linspace(1,30,30)*1e-2
FF=np.linspace(1,10,2)*1e-2
delta=np.zeros(FF.shape[0])
R0=np.zeros(FF.shape[0])
RI=np.zeros(FF.shape[0])
R1=np.zeros(FF.shape[0])
for i,F_inp in enumerate(FF):
   R0[i], RI[i], R1[i], delta[i] = force_dist(F_inp,t_R0,t_R1)
   t_R0=np.arange(R0[i],R0[i]+0.05,0.01)
   t_R1=np.arange(R1[i], R1[i]+0.02, 1e-3)
dd=np.vstack([delta,FF,R0,RI,R1])
np.savetxt("data/sigma0o1_k0o1.txt",dd)
#+end_src
#+RESULTS: