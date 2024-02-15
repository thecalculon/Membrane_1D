import sys
from Membrane_1D import *
import matplotlib.pyplot as plt
#--------------------------------------------------------------#
def plot_mem(ax,fname,point):
    Np=4096
    dd=np.loadtxt(fname)
    rr,zz,rrtip,zztip=membrane(dd[point,2],dd[point,3],dd[point,4],Np=Np)
    print((dd[point,0]-dd[0,0])*radius_ev*nano,
        (dd[point,1]-dd[0,1])*radius_ev*pico)
    plt.plot(rr[0:],zz[0:],'k-',linewidth=1.0)
    plt.plot(-rr[4*Np-1::-1],zz[4*Np-1::-1],'k-',linewidth=1.0)
    plt.plot(rrtip,zztip,'r-',linewidth=1.0)
    plt.plot(-rrtip[-1::-1],zztip[-1::-1],'r-',linewidth=1.0)
    ax.set_aspect('equal')
    ax.spines[['right','top','left','bottom']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
#--------------------------------------------------------------#
fname=sys.argv[1]
#--------------------------------------------------------------#
fig,ax=plt.subplots(1,1)
plot_mem(ax,fname+".txt",1)
plot_mem(ax,fname+".txt",25)
plot_mem(ax,fname+".txt",80)
# plt.savefig("../figures/"+fname.replace("data","")+"_geomtry.pdf")
plt.show()