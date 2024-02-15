# Make the membrane.
import numpy as np
#################################################################
Rv = 1;
Av = 4*np.pi*Rv*Rv;
VolT = 4./3.*np.pi* Rv*Rv*Rv
#################################################################
def make_membrane(Np=128, radius=1.0):
    dx = 2*np.pi/(Np)
    xx = np.linspace(0,2*np.pi,Np)
    xp = radius*np.sin(xx) 
    yy = radius*np.cos(xx)
    return xp, yp
#--------------------------------------------------------------#
def a1b1(R0, RI, R1):
    A1 = R0/(R0*R0 - RI*RI)
    B1 = -A1*RI*RI
    return A1, B1
#--------------------------------------------------------------#
def a3b3(R0, RI, R1, theta=0.26):
    A3 = R1*np.cos(theta)/(RI*RI - R1*R1)
    B3 = -A3*R1*R1 - R1*np.cos(theta)
    return A3, B3
#--------------------------------------------------------------#
def __u(A, B, r):
    if(r<1e-5):
        print("Exception raise: r=0")
    else:
        u = A*r + B/r
    return u
#--------------------------------------------------------------#
def Area(R0, RI, R1, Np=8192, theta=0.26):
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
#--------------------------------------------------------------#
def Volume(R0, RI, R1, Np=4096, theta=0.26):
    dr = (R0 - RI)/Np
    rr = np.arange(RI+dr/2, R0+dr/2, dr)
    A1, B1 = a1b1(R0, RI, R1)
    A3, B3 = a3b3(R0, RI, R1)
    vol1, vol2, vol3 = 0.0, 0.0, 0.0;
    t1, t2, t3 = 0.0, 0.0, 0.0;
    for r in rr:
        u1 = __u(A1, B1, r)
        den1 = 1.0 - u1*u1
        if (den1 > 1e-12):
           t1 += u1*r*r*dr/np.sqrt(den1);
           t3 += u1*dr/np.sqrt(den1)
    vol1 = np.pi*(t1 + t1 - RI*RI*t3)
    rr = np.arange(R1+dr/2, RI+dr/2, dr)
    for r in rr:
        u3 = __u(A3, B3, r)
        den2 = 1.0 - u3*u3
        if (den2 > 1e-12):
           t1 = u3*r*r*dr/np.sqrt(den2);
           vol2 += np.pi*t1
    vol3 = -np.pi*R1*R1*R1/(3*np.tan(theta))
    return vol1, vol2, vol3
#--------------------------------------------------------------#
def force(R0, RI, R1, kasigma, Np=128, Aind=None, theta=0.26):
   if Aind is None:
      At, Ab = Area(R0, RI, R1, Np=128)
      Aind = At + Ab
   A3, B3 = a3b3(R0, RI, R1)
   t1 = 2*np.pi*(R1*np.sin(0.5*np.pi - theta) + R1*R1*A3)
   t2 = kasigma[1]+kasigma[0]*(Aind - Av)/Av
   ft = t1*t2
   fb = ft
   return ft, fb
#--------------------------------------------------------------#
def delta(R0, RI, R1, Np=128, theta=0.26):
    dr = (R0 - RI)/Np
    A1, B1 = a1b1(R0, RI, R1)
    A3, B3 = a3b3(R0, RI, R1)
    t1, t2 = 0.0, 0.0;
    rr = np.arange(RI+dr/2, R0+dr/2, dr)
    for r in rr:
        u1 = __u(A1, B1, r)
        denom = 1-u1*u1
        if (denom > 1e-12):
            t1+=u1*dr/np.sqrt(denom)
    #
    rr = np.arange(R1+dr/2, RI+dr/2, dr)
    for r in rr:
        u3 = __u(A3, B3, r)
        denom=1-u3*u3
        if (denom > 1e-12):
            t2+=u3*dr/np.sqrt(denom)
    delta= 2*Rv - (2*t1 + t2 - R1/np.tan(theta))
    return delta
#--------------------------------------------------------------#
def getRI(R0,R1,theta=0.26):
    num=R0*R0*R1*np.cos(theta)+R0*R1*R1
    denom=R0+R1*np.cos(theta)
    RI = np.sqrt(num/denom)
    return RI
#--------------------------------------------------------------#
def tip(Rst,Np=128,theta=0.26):
    dr=Rst/Np
    rr=np.arange(Rst,0,-dr)
    zz=rr/np.sin(theta)
    return rr,zz
#--------------------------------------------------------------#
def forcedistcurve(kasigma,delta_max=None,Np=100,verbose=True,
    theta=0.26):
    t_R1=np.linspace(1e-6,12e-2,Np)
    t_R0=np.arange(1,1.02,0.0001)
    t_delta=np.empty(Np,dtype='float64')
    t_FF=np.empty(Np,dtype='float64')
    i=0
    t_R0_fin=np.zeros(Np)
    t_RI=np.zeros(Np)
    for R1 in t_R1:
        for R0 in t_R0:
            RI=getRI(R0,R1)
            V1, V2, V3 = Volume(R0, RI, R1, Np=4096)
            Vol = V1+V2+V3
            At, Ab=Area(R0, RI, R1)
            cdt4 = At+Ab>Av
            cdt3 = abs(Vol - VolT) < (1e-3)*VolT;
            if cdt3 and cdt4:
                t_delta[i]=delta(R0,RI,R1)
                t_FF[i]=force(R0, RI, R1, kasigma, Aind=At+Ab)[0]
                t_R0_fin[i]=R0
                t_RI[i]=RI
                i=i+1
                break;
        if verbose:
            print(t_delta[i-1],t_FF[i-1])
        if(delta_max is not None and t_delta[i-1]-t_delta[0]>delta_max):
            break;
        t_R0=np.arange(R0,R0+0.001,0.0001)
    return t_delta[0:i],t_FF[0:i],t_R0_fin[0:i],t_RI[0:i],t_R1[0:i]
#--------------------------------------------------------------#
def chop2(xed,meany,min_=None,max_=None):
    if min_ is None:
        min_=0
    ind=xed>min_
    xed,meany=xed[ind],meany[ind]
    if max_ is not None:
        ind=xed<max_
        xed,meany=xed[ind],meany[ind]
    return xed,meany
#--------------------------------------------------------#
def chop(xed,meany,meanperr,meanmerr,min_=None,max_=None):
    if min_ is None:
        min_=0
    ind=xed>min_
    xed,meany,meanperr,meanmerr=xed[ind],meany[ind],meanperr[ind],meanmerr[ind]
    if max_ is not None:
        ind=xed<max_
        xed,meany,meanperr,meanmerr=xed[ind],meany[ind],meanperr[ind],meanmerr[ind]
    return xed,meany,meanperr,meanmerr
#--------------------------------------------------------#