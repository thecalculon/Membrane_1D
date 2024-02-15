import numpy as np
# from func import *
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
def membrane(R0,RI,R1,Np=128,theta=0.26):
    zz=np.zeros(Np*4)
    rr=np.zeros(Np*4)
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
    dr=R1/Np
    rrtip1=np.arange(rr[4*Np-1],0,-dr)
    zztip1=zz[4*Np-1]-(rr[4*Np-1]-rrtip1)/np.tan(theta)
    rrtip2=np.arange(rr[4*Np-1],2*rr[4*Np-1],dr)
    zztip2=zz[4*Np-1]+(rrtip2-rr[4*Np-1])/np.tan(theta)
    return rr,zz,np.hstack([rrtip1,rrtip2]),np.hstack([zztip1,zztip2])
#--------------------------------------------------------------#