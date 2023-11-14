import numpy as np
import matplotlib.pyplot as plt
import sys
# dd=np.loadtxt("data/sigma0o001_k1.txt")
dd=np.loadtxt(sys.argv[1])
plt.plot(dd[0,1:],dd[1,1:],'o-')
# dd=np.loadtxt(sys.argv[2])
# plt.plot(dd[0,1:],dd[1,1:],'*-')
plt.show()