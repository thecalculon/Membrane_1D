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
#+begin_src python :session mempy :results output
#+end_src
#+begin_src python :session mempy :results output
#################################################################
# Lets start playing around
# R0>RI>R1
F_inp = 0.3
t_R0=np.arange(0.8, 2.0, 0.01)
t_RI=np.arange(0.02, 0.79, 0.01)
t_R1=np.arange(0.002, 0.019, 0.01)
for R0 in t_R0:
   for RI in t_RI:
      for R1 in t_R1:
        print("Checking for R0=", R0, ", RI=", RI, ", R1=", R1)
        ft, fb =  force_all(R0, RI, R1, Np=4096)
        V1, V2, V3 =  Volume(R0, RI, R1, Np=4096)
        Vol = V1+V2+V3
        cdt1 = abs(ft - F_inp) < 0.05;
        cdt2 = abs(fb - F_inp) < 0.05;
        cdt3 = abs(Vol - VolT) < 0.1;
        if(cdt1 and cdt2 and cdt3): print(R0, RI, R1)
        if(cdt1 and cdt2): print("1 and 2", ft, fb, Vol, R0, RI, R1)
        if(cdt2 and cdt3): print("2 and 3", ft, fb, Vol, R0, RI, R1)
        if(cdt1 and cdt3): print("1 and 3", ft, fb, Vol, R0, RI, R1)
#+end_src
#+RESULTS: