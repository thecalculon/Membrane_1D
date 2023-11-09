from func import *
import matplotlib.pyplot as plt
#--------------------------------------------------------------#
def membrane2(R0,RI,R1,Np=128):
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
# dd=np.loadtxt("data/sigma0o1_k0o1.txt")
# rr,zz=membrane(dd[2,10],dd[3,10],dd[4,10])
lambd=0.2
R0=1
rr,zz=membrane(R0,lambd,0)
Np=int(zz.shape[0]/8)
fig, ax = plt.subplots(1, 1)
plt.plot(rr[0:],zz[0:],'o')
ax.set_aspect('equal')
plt.show()