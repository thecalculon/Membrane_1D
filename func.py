# Make the membrane.
from header import *
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
def a2b2(R0, RI, R1):
    A2, B2 = a1b1(R0, RI, R1)
    return A2, B2
#--------------------------------------------------------------#
def a3b3(R0, RI, R1):
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
def Area(R0, RI, R1, Np=8192):
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
#--------------------------------------------------------------#
def Volume(R0, RI, R1, Np=4096):
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
        # if (den2 > 1e-12):
           # t2 += u3*r*r*dr/np.sqrt(den2);
    vol1 = np.pi*(t1 + t1 - RI*RI*t3)
    # print(t1, t2, RI*RI*t3)
    rr = np.arange(R1+dr/2, RI+dr/2, dr)
    for r in rr:
        u3 = __u(A3, B3, r)
        den2 = 1.0 - u3*u3
        if (den2 > 1e-12):
           t1 = u3*r*r*dr/np.sqrt(den2);
           vol2 += np.pi*t1
    vol3 = - np.pi*R1*R1*R1/(3*np.tan(theta))
    return vol1, vol2, vol3
#--------------------------------------------------------------#
def force(R0, RI, R1, kasigma, Np=128, Aind=None):
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
def delta(R0, RI, R1, Np=128):
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
def getRI(R0,R1):
    num=R0*R0*R1*np.cos(theta)+R0*R1*R1
    denom=R0+R1*np.cos(theta)
    RI = np.sqrt(num/denom)
    return RI
#--------------------------------------------------------------#
def getzz(rr,z0,A,B):
    t1=z0
    dr=rr[1]-rr[0]
    zz=np.zeros(rr.shape[0])
    for ip,r in enumerate(rr):
        u=__u(A,B,r)
        t1+=u*dr/np.sqrt(1-u*u)
        zz[ip]=t1
    return zz
#--------------------------------------------------------------#
def membrane(R0,RI,R1,Np=128):
    zz=np.zeros(Np*8)
    rr=np.zeros(Np*8)
    A1,B1=a1b1(R0, RI, R1)
    A2,B2=A1,B1
    A3,B3=a3b3(R0, RI, R1)
    A3,B3=-A3,-B3
    ##
    dr=(RI-0)/Np
    rr[0:Np]=np.arange(0,RI,dr)
    zz[0:Np]=0
    ##
    dr=(R0-RI)/Np
    rr[Np:2*Np]=np.arange(RI,R0,dr)
    zz[Np:2*Np]=getzz(rr[Np:2*Np],zz[Np-1],A1,B1)
    ##
    rr[2*Np:3*Np]=rr[2*Np-1:Np-1:-1]
    zz[2*Np:3*Np]=2*zz[2*Np-1]-zz[2*Np-1:Np-1:-1]
    ##
    dr=(RI-R1)/Np
    rr[3*Np:4*Np]=np.arange(RI,R1,-dr)
    zz[3*Np:4*Np]=getzz(rr[3*Np:4*Np],zz[3*Np-1],A3,B3)
    ##
    rr[4*Np:]=-rr[4*Np-1::-1]
    zz[4*Np:]=zz[4*Np-1::-1]
    ##
    return rr,zz
#--------------------------------------------------------------#
def delta4force(F_inp,kasigma,t_R0,t_R1):
    for R0 in t_R0:
        for R1 in t_R1:
            RI=getRI(R0,R1)
            ft, fb = force(R0, RI, R1, kasigma, Np=4096)
            V1, V2, V3 =  Volume(R0, RI, R1, Np=4096)
            Vol = V1+V2+V3
            # print(ft, fb)
            print(R0,R1,(ft-F_inp)/F_inp, (Vol - VolT)/VolT)
            cdt1 = abs(ft - F_inp) < 0.05*F_inp;
            cdt3 = abs(Vol - VolT) < 0.05*VolT;
            if(cdt1 and cdt3):
                print("ft=", ft, "Vol=", Vol, "R0=", R0, "RI=", RI, "R1=", R1)
                return R0, RI, R1, delta(R0, RI, R1)
#--------------------------------------------------------------#
def forcedistcurve(kasigma,delta_max=None,Np=100,verbose=True):
    t_R1=np.linspace(1e-6,5e-2,Np)
    t_R0=np.arange(1,1.02,0.0001)
    t_delta=np.empty(Np,dtype=object)
    t_FF=np.empty(Np,dtype=object)
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
            print(t_delta[i-1],t_FF[i-1],R0,RI,R1)
        if(delta_max is not None and t_delta[i-1]-t_delta[0]>delta_max):
            break;
        t_R0=np.arange(R0,R0+0.001,0.0001)
    t_delta[0:i]=t_delta[0:i]-t_delta[0]
    return t_delta[0:i],t_FF[0:i],t_R0_fin[0:i],t_RI[0:i],t_R1[0:i]
#--------------------------------------------------------------#
def cost_func(kasigma,expt_data,Np=100):
    kasigma=kasigma.T
    print(kasigma)
    t_delta,t_FF=forcedistcurve(kasigma,np.max(expt_data[:,0]),Np,False)
    t_delta=t_delta-t_delta[0]
    area1=np.trapz(t_FF,t_delta)
    area2=np.trapz(expt_data[:,1],expt_data[:,0])
    print("cost=",np.abs(area2-area1))
    cost=np.abs(area2-area1)/(kasigma[0]*Rv*Rv)
    return cost
#--------------------------------------------------------------#