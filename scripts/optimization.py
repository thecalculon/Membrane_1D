import numpy as np
from Membrane_1D import *
import pandas as pd
import numpy.linalg as la
import glob
from scipy.optimize import curve_fit
import os
#################################################################
frcecomp=False
onebynano=1e9
onebypico=1e12
theta=0.26
radius_ev=84.38029625e-9
#################################################################
def cost(xx,ka,sigma):
    sigma=round(sigma,11)
    ka=round(ka,9)
    kasigma=np.array([ka,sigma])
    fname="data/sigma"+str(sigma).replace(".","o")+\
    "_k"+str(ka).replace(".","o")+"_th"+\
    str(theta).replace(".","o")+".txt"
    if os.path.isfile(fname) and frcecomp==False:
        print("Taking data from "+ fname)
        dd=np.loadtxt(fname)
        t_delta,t_FF=dd[:,0],dd[:,1]
    else:
        print("generating data for Ka =", str(ka), ", sigma =", str(sigma))
        t_delta,t_FF,a,b,c=forcedistcurve(kasigma,np.max(xx)/(radius_ev*onebynano),
			verbose=False)
        np.savetxt(fname,np.vstack([t_delta,t_FF,a,b,c]).T)
    t_delta=(t_delta-t_delta[0])*radius_ev*onebynano
    t_FF=t_FF*radius_ev*onebypico
    FF_model=np.interp(xx,t_delta,t_FF)
    return np.log10(FF_model)
#################################################################
vesicle_fols=['wt','cd63ko','panko','wt_sec']
for j in range(1):
    dirs=sorted(glob.glob("/home/vipin/PKO_expt/AFM_EV/"+vesicle_fols[j]+"/*"))
    vesdata=pd.read_excel("/home/vipin/PKO_expt/AFM_EV/vesicle_data_other_RC.xlsx",
                        sheet_name=j)
    df1 = pd.DataFrame([['Vesicle', 'Ka', 'Signot']], 
    columns=['Vesicle', 'Ka', 'Signot'])
    for i,fol in enumerate(dirs):
        vesname=fol.split("/")[-1]
        xmax=list(vesdata['Xmax'][vesdata['Name']==vesname])[0]
        radius_ev=list(vesdata['Rc_kbt'][vesdata['Name']==vesname])[0]*1e-9
        if xmax>0:
            print("\nOptimizing the data in ",fol)
            dd=np.loadtxt(fol+"/mean_fd.txt")
            dexpt,fexpt=dd[:,0]*onebynano,dd[:,1]*onebypico
            dexpt,fexpt=chop2(dexpt,fexpt,min_=0,max_=xmax)
            dexpt=dexpt+list(vesdata['Shift'][vesdata['Name']==vesname])[0]
            d1,fpush=dexpt,np.log10(fexpt)
            popt, pcov = curve_fit(cost, d1, fpush, p0=[0.5,0.005],
                method='lm', xtol=1e-2)
            print("popt=",popt)
            print("Error =", np.sqrt(np.dot(10**fpush-10**cost(d1,*popt),
                10**fpush-10**cost(d1,*popt))/np.dot(10**fpush,10**fpush))*100)
            print("condition = ", la.cond(pcov))
            df2=pd.DataFrame([[vesname,popt[0],popt[1]]],
                columns=['Vesicle', 'Ka', 'Signot'])
            df1=pd.concat([df1,df2],ignore_index=True)
    df1.to_csv(r'Ka_'+vesicle_fols[j], header=None, index=None, sep=' ')
#-------------------------------------------------------------#