# Following are the assumptions
# R1<RI<R0
# R0>RI>R1
# The only issue is: r = 0 seems singular

# Coordinate system::
# ^ z
# |
# ------> r


#+begin_src python :session mempy :results output

import numpy as np
import matplotlib.pyplot as plt

#+end_src


#+begin_src python :session mempy :results output

# Make the membrane.
def make_membrane(Np=128, radius=1.0):
    dx = 2*np.pi/(Np)
    xx = np.linspace(0,2*np.pi,Np)
    xp = radius*np.sin(xx) 
    yy = radius*np.cos(xx)
    return xp, yp

def a1b1(R0, RI, R1):
    A1 = R0/(R0*R0 - RI*RI)
    B1 = -A1*RI*RI
    return A1, B1

def a2b2(R0, RI, R1):
    A2, B2 = a1b1(R0, RI, R1)
    return A2, B2

def a3b3(R0, RI, R1):
    A3 = R1*np.sin(0.5*np.pi - theta)/(RI*RI - R1*R1)
    B3 = -A3*R1*R1 - R1*np.sin(0.5*np.pi - theta)
    return A3, B3

def __u(A, B, r):
    if(r<1e-5):
       print("Exception raise: r=0")
    else:
       u = A*r + B/r
    return u

#+end_src




#+begin_src python :session mempy :results output

# Here are the functions for estimation of area


def Area(R0, RI, R1, Np=128):
    dr = (R0 - RI)/Np
    rr = np.arange(RI+dr/2, R0+dr/2, dr)
    A1, B1 = a1b1(R0, RI, R1)
    A3, B3 = a3b3(R0, RI, R1)
    Arbot1, Artop1, Artop2 = 0.0, 0.0, 0.0;
    for r in rr:
        u = __u(A1, B1, r)
        denom = 1.0 - u*u
        if (denom > 1e-12):
           Arbot1 += r*dr/np.sqrt(denom);
        # Artop1 += r*dr/np.sqrt(1 - u*u);
    Artop1 = Arbot1;

    rr = np.arange(R1+dr/2, RI+dr/2, dr)
    for r in rr:
        u = __u(A3, B3, r)
        denom = 1.0 - u*u
        if (denom > 1e-12):
           Artop2 += r*dr/np.sqrt(denom);
 
    Arbot = np.pi*(RI*RI + 2*Arbot1)
    Artop = np.pi*(2*Artop1 + 2*Artop2 + R1*R1/np.sin(theta))

    return Arbot, Artop


def Volume(R0, RI, R1, Np=128):
    dr = (R0 - RI)/Np
    rr = np.arange(RI+dr/2, R0+dr/2, dr)
    A1, B1 = a1b1(R0, RI, R1)
    A3, B3 = a3b3(R0, RI, R1)
    vol1, vol2, vol3 = 0.0, 0.0, 0.0;
    t1, t2, t3 = 0.0, 0.0, 0.0;
    for r in rr:
        u1 = __u(A1, B1, r)
        u3 = __u(A3, B3, r)
        den1 = 1.0 - u1*u1
        den2 = 1.0 - u3*u3
        if (den1 > 1e-12):
           t1 += u1*r*r*dr/np.sqrt(den1);
           t3 += u1*dr/np.sqrt(den1) 
        if (den2 > 1e-12):
           t2 += u3*r*r*dr/np.sqrt(den2);

    vol1 = np.pi*(t1 + t2 - RI*RI*t3) 
    # print(t1, t2, RI*RI*t3)
    rr = np.arange(R1+dr/2, RI+dr/2, dr)
    for r in rr:
        u3 = __u(A3, B3, r)
        den2 = 1.0 - u3*u3
        if (den2 > 1e-12):
           t1 = u3*r*r*dr/np.sqrt(den2);
           vol2 += np.pi*t1

    vol3 = - np.pi*R1*R1*R1/(3*np.tan(theta))
    vol = vol1 + vol2 + vol3 # - np.pi*R1*R1*R1/(3*np.tan(theta))

    return vol1, vol2, vol3



#+end_src





#+begin_src python :session mempy :results output

# here is the function for force balance
def force_all(R0, RI, R1, Np=128):
    At, Ab = Area(R0, RI, R1, Np=128)
    Aind = At + Ab
    A1, B1 = a1b1(R0, RI, R1)
    A3, B3 = a3b3(R0, RI, R1)
    t1 = 2*np.pi*(R1*np.sin(0.5*np.pi - theta) + R1*R1*A3)
    t2 = sigma_0 + KA*(Aind - Av)/Av
    ft = t1*t2
    fb = t2*2*np.pi*RI*RI*A1 
    return ft, fb


#+end_src




#+begin_src python :session mempy :results output

# Lets start playing around
sigma_0 = 0.1; KA = 0.1; theta = 0.31; 
Rv = 1.0; Av = 4*np.pi*Rv*Rv; VolT = 4./3. * np.pi* Rv*Rv*Rv
# R0>RI>R1
F_inp = 0.3
t_R0 = np.arange(0.8, 2.0, 0.01 )
t_RI = np.arange(0.02, 0.7, 0.01 )
t_R1 = np.arange(0.002, 0.019, 0.01 )
for R0 in t_R0:
   for RI in t_RI:
      for R1 in t_R1:
          ft, fb =  force_all(R0, RI, R1, Np=1024)
          V1, V2, V3 =  Volume(R0, RI, R1, Np=1024)
          Vol = V1+V2+V3
          cdt1 = abs(ft - F_inp) < 0.05;
          cdt2 = abs(fb - F_inp) < 0.05;
          cdt3 = abs(2*Vol - VolT) < 0.1;
          if(cdt1 and cdt2 and cdt3): print(R0, RI, R1)
          if(cdt1 and cdt2): print("1 and 2", ft, fb, Vol, R0, RI, R1)
          if(cdt2 and cdt3): print("2 and 3", ft, fb, Vol, R0, RI, R1)
          if(cdt1 and cdt3): print("1 and 3", ft, fb, Vol, R0, RI, R1)

#+end_src

#+RESULTS:

