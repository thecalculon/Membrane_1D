import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/home/vipin/.matplotlibrc')
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
#######################################################
sigma_0 = 1e-4;
KA = 5e-1;
theta = 0.314;
Rv = 1;
Av = 4*np.pi*Rv*Rv;
VolT = 4./3.*np.pi* Rv*Rv*Rv
fname="data/sigma"+str(sigma_0).replace(".","o")+\
	"_k"+str(KA).replace(".","o")+".txt"
print(fname)
#######################################################