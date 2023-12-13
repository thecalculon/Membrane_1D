import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/home/vipin/.matplotlibrc')
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
import lib as lib
import glob as glob
import os
#######################################################
sigma_0 = 1
KA = 1;
theta = 0.18;
Rv = 1;
Av = 4*np.pi*Rv*Rv;
VolT = 4./3.*np.pi* Rv*Rv*Rv
fname="data/sigma"+str(sigma_0).replace(".","o")+\
	"_k"+str(KA).replace(".","o")+".txt"
compute=False
print(fname)
#######################################################
kasigma=np.array([KA,sigma_0])
expt_data=np.loadtxt("data/sigma0o0001_k0o5.txt").T
expt_data=expt_data[0:20,:]
#-----------------------------------------------------#
kasigmax=np.array([1,1])
kasigmin=np.array([0,0])
#######################################################