import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import pandas as pd
import lib as lib
import matplotlib
import matplotlib.cm as clr_b
######################################################################
nano=1e9
######################################################################
def sort_(fnames):
    tidx = [f[-16:-4].replace('.', '') for f in fnames]
    tidx = np.asarray(tidx, dtype=int)
    idx = np.argsort(tidx)
    fn = [fnames[i] for i in idx]
    return fn
#--------------------------------------------------------------------#
def pickled_data(fol):
    files_ = glob.glob(fol+"/force-save*")
    files = sort_(files_)
    print(files)
    df_push, tmp_df = pd.DataFrame({}), pd.DataFrame({})
    df_pull, tmp_df = pd.DataFrame({}), pd.DataFrame({})
    d1, d2, fpush, fretrace = lib.read_file(files[1]);
    char = "%04d" %(0)
    df_push[char] = np.vstack([d1, fpush]).tolist()
    df_pull[char] = np.vstack([d2, fretrace]).tolist()
    #
    for i, f in enumerate(files[2:]):
        char = "%04d" %(i+1)
        dt1, dt2, ftpush, ftretrace = lib.read_file(f);
        tmp_df = pd.DataFrame()
        tmp_df[char] = np.vstack([dt1, ftpush]).tolist()
        df_push = pd.concat([df_push, tmp_df], axis=1)
        tmp_df = pd.DataFrame()
        tmp_df[char] = np.vstack([dt2, ftretrace]).tolist()
        df_pull = pd.concat([df_pull, tmp_df], axis=1)
    df_push.to_pickle(fol+'/force_push.pkl')
    df_pull.to_pickle(fol+'/force_pull.pkl')
#--------------------------------------------------------------------#
def moving_average(dataset,window_size=5):
	result = []
	for i in range(len(dataset) - window_size + 1):
	    window = dataset[i : i + window_size]
	    window_average = sum(window) / window_size
	    result.append(window_average)
	return np.array(result)
#--------------------------------------------------------------------#
def mean_indexp(fol,df_push=None,startiter=0):
    if df_push is None:
        df_push = pd.read_pickle(fol+"/force_push.pkl")
    df_pull = pd.read_pickle(fol+"/force_pull.pkl")
    dpush, dpull, fpush, fpull = lib.stacked_force_dis(df_push, df_pull, df_push.keys())
    nbins = 128
    H, xedges, yedges = np.histogram2d(dpush, fpush, bins=nbins, density = True)
    xed = 0.5*(xedges[0:-1] + xedges[1:])
    yed = 0.5*(yedges[0:-1] + yedges[1:])
    H = H.T
    xedn, yedn = np.meshgrid(xed, yed)
    idx = H.argmax(axis=0)
    meany = np.zeros(len(yed))
    sigma = np.zeros(len(yed))
    for i in range(len(yed)):
        hist = H[:,i]/(np.sum(H[:,i]))
        meany[i] = np.sum(hist*yed)
        sigma[i] = np.sum(hist*yed*yed)
    #
    sigma = sigma - meany*meany
    meanperr = meany+np.sqrt(sigma)
    meanmerr = meany-np.sqrt(sigma)
    return xed,meany,meanperr,meanmerr
#--------------------------------------------------------------------#
def density_plot_iterator(fol,startiter=0,savetxt=False):
    df_push = pd.read_pickle(fol+"/force_push.pkl")
    xed,meany,meanperr,meanmerr=mean_indexp(fol,df_push,startiter=startiter)
    fig,ax=plt.subplots()
    time_plot(fig, ax, df_push)
    ax.plot(xed*nano, meany*nano, '-', color='tab:red', linewidth=2)
    ax.plot(xed*nano, meanperr*nano, '--k', linewidth=2)
    ax.plot(xed*nano, meanmerr*nano, '--k', linewidth=2)
    ax.grid(color='k', linestyle='--', linewidth=0.5)
    if savetxt==True:
        np.savetxt(fol+'/mean_fd.txt', np.vstack([xed,meany,meanperr]).T)
#--------------------------------------------------------------------#
def time_plot(fig, ax, push):
    numkey = len(push.keys())
    norm = matplotlib.colors.Normalize(vmin = 0,vmax = numkey)
    for i, a in enumerate(push.keys()):
        rgba = clr_b.viridis(norm(i), bytes=False)
        rgba_hex = matplotlib.colors.rgb2hex(rgba)
        d = np.asarray(push[a][0])*nano
        force = np.asarray(push[a][1])*nano
        ax.plot(-d, force, '-', color=rgba_hex)
#--------------------------------------------------------------------#
def chop(xed,meany,meanperr,meanmerr,min_=None,max_=None):
    if min_ is None:
        min_=0
    ind=xed>min_    
    xed,meany,meanperr,meanmerr=xed[ind],meany[ind],meanperr[ind],meanmerr[ind]
    if max_ is not None:
        ind=xed<max_
        xed,meany,meanperr,meanmerr=xed[ind],meany[ind],meanperr[ind],meanmerr[ind]
    return xed,meany,meanperr,meanmerr
#--------------------------------------------------------------------#
def chop2(xed,meany,min_=None,max_=None):
    if min_ is None:
        min_=0
    ind=xed>min_    
    xed,meany=xed[ind],meany[ind]
    if max_ is not None:
        ind=xed<max_
        xed,meany=xed[ind],meany[ind]
    return xed,meany
######################################################################
fol="Ref_data/WT/221222ev5_h_70nm"
# pickled_data("Ref_data/WT/221222ev5_h_70nm/")
xed,meany,meanperr,meanmerr=mean_indexp(fol)
# xed,meany,meanperr,meanmerr=chop(xed,meany,meanperr,meanmerr,min_=0,max_=2.2e-8)
np.savetxt(fol+"/mean_fd.txt",np.vstack([xed,meany]).T)
fig,ax=plt.subplots()
ax.plot(xed*nano, meany*nano, 'o-', color='tab:red', linewidth=2)
#
# files=glob.glob("Ref_data/fdc_example/force-save*")
# d1,d2,fpush,fret=lib.read_file(files[53])
# print(files[53])
# d1,fpush=chop2(-d1,fpush,min_=-2e-9,max_=2.2e-8)
# d1_sm=moving_average(d1,75)
# fpush_sm=moving_average(fpush,75)
# #
# ax.loglog(d1_sm*nano,fpush_sm*nano,'.', label='window=50')
plt.show()