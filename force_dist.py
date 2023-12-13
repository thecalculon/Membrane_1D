from header import *
from func import *
#########################################################
nano=1e9
radius_ev=80.4e-9
#########################################################
def moving_average(dataset,window_size=5):
	result = []
	for i in range(len(dataset) - window_size + 1):
	    window = dataset[i : i + window_size]
	    window_average = sum(window) / window_size
	    result.append(window_average)
	return np.array(result)
#--------------------------------------------------------#
def chop(xed,meany,min_=None,max_=None):
    if min_ is None:
        min_=0
    ind=xed>min_    
    xed,meany=xed[ind],meany[ind]
    if max_ is not None:
        ind=xed<max_
        xed,meany=xed[ind],meany[ind]
    return xed,meany
#########################################################
fig,ax=plt.subplots()
# files=glob.glob("Ref_data/fdc_example/force-save*")    
# d1,d2,fpush,fret=lib.read_file(files[53])
# d1,fpush=chop(d1,fpush,min_=0,max_=2e-8)
# d1_sm=moving_average(d1,75)
# fpush_sm=moving_average(fpush,75)
dd=np.loadtxt("Ref_data/WT/0616ev402.txt")
# np.savetxt("Ref_data/GVB/fig12p5a_sigma0o0001_k0o25.csv",dd)
d1,fpush=dd[:,0],dd[:,1]
d1,fpush=chop(d1,fpush,min_=0)
# d1_sm,fpush_sm=chop(d1,fpush,min_=0,max_=2e-8)
#----------------------------------------------------#
if os.path.isfile(fname) and compute==False:
    dd=np.loadtxt(fname)
    t_delta,t_FF=dd[:,0],dd[:,1]
else:
    t_delta,t_FF,t_R0,t_RI,t_R1=forcedistcurve(kasigma,delta_max=0.35)
    np.savetxt(fname,np.vstack([t_delta,t_FF,t_R0,t_RI,t_R1]).T)
#----------------------------------------------------#
# d1=d1-0.058*radius_ev
# fpush_sm=fpush_sm
ax.plot(d1/radius_ev,fpush/(radius_ev*KA),'o-')
ax.plot(t_delta[1:],t_FF[1:]/KA,'o-')
plt.show()