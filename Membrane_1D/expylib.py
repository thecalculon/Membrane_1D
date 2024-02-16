import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import glob as glob
import pandas as pd

def sort_(fnames):
    tidx = [f[-16:-4].replace('.', '') for f in fnames]
    # print(fnames,tidx)
    tidx = np.asarray(tidx, dtype=int)
    idx = np.argsort(tidx)
    fn = [fnames[i] for i in idx]
    return fn

def read_file(f):
    """
    f: filename
    returns: data read from the ascii file
    """
    d1, d2, fpush, fretrace = np.asarray([]), np.asarray([]), np.asarray([]), np.asarray([])
    f = open(f, "r")
    lines = f.readlines()
    is_push = True
    for line in lines: #.split("\n"):
        if line == '\n':
            is_push = False
        if is_push and not line.startswith('#'):
            d = line.strip().split()
            d1 = np.hstack([d1,float(d[0])])
            fpush = np.hstack([fpush,float(d[1])])
        if line != '\n' and not is_push and not line.startswith('#'):
            d = line.strip().split()
            d2 = np.hstack([d2,float(d[0])])
            fretrace = np.hstack([fretrace,float(d[1])])
    return d1, d2, fpush, fretrace


def trunc_fd(push, pull):
    """
    push pull: are the pd data for certain keys
    pull curve hovers over the starting point a lot. It truncates the pull curve
    returns: force-push, force-pull, push_d, pull_d
    """
    fpush, fpull = np.asarray(push[1]), np.asarray(pull[1])
    disps, displ = np.asarray(push[0]), np.asarray(pull[0])
    fpull = fpull[0:len(fpush)]
    displ = displ[0:len(fpush)]
    return fpush, fpull, -disps, -displ

def stacked_force_dis(push, pull, keys):
    dpush, dpull = np.asarray([]), np.asarray([])
    fpush, fpull = np.asarray([]), np.asarray([])
    for a in keys:
        tfpush, tfpull, tdpush, tdpull = trunc_fd(push[a], pull[a])
        dpush = np.hstack([dpush, tdpush])
        dpull = np.hstack([dpull, tdpull])
        fpush = np.hstack([fpush, tfpush])
        fpull = np.hstack([fpull, tfpull])
    return dpush, dpull, fpush, fpull

def pickled_data(fol):
    files_ = glob.glob(fol+"/after-fdc*")
    print(files_)
    files = sort_(files_)
    df_push, tmp_df = pd.DataFrame({}), pd.DataFrame({})
    df_pull, tmp_df = pd.DataFrame({}), pd.DataFrame({})
    d1, d2, fpush, fretrace = read_file(files[1]);
    char = "%04d" %(0)
    df_push[char] = np.vstack([d1, fpush]).tolist()
    df_pull[char] = np.vstack([d2, fretrace]).tolist()
    #
    for i, f in enumerate(files[2:]):
        char = "%04d" %(i+1)
        dt1, dt2, ftpush, ftretrace = read_file(f);
        tmp_df = pd.DataFrame()
        tmp_df[char] = np.vstack([dt1, ftpush]).tolist()
        df_push = pd.concat([df_push, tmp_df], axis=1)
        tmp_df = pd.DataFrame()
        tmp_df[char] = np.vstack([dt2, ftretrace]).tolist()
        df_pull = pd.concat([df_pull, tmp_df], axis=1)
    df_push.to_pickle(fol+'/force_push.pkl')
    df_pull.to_pickle(fol+'/force_pull.pkl')

def mean_indexp(fol,df_push=None,startiter=0):
    if df_push is None:
        df_push = pd.read_pickle(fol+"/force_push.pkl")
    df_pull = pd.read_pickle(fol+"/force_pull.pkl")
    dpush, dpull, fpush, fpull = lib.stacked_force_dis(df_push, df_pull, df_push.keys())
    nbins = 128
    H, xedges, yedges = np.histogram2d(dpush[dpush>-2e-7], fpush[dpush>-2e-7], bins=nbins, density = True)
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
    sigma = sigma - meany*meany
    meanperr = meany+np.sqrt(sigma)
    meanmerr = meany-np.sqrt(sigma)
    return xed,meany,meanperr,meanmerr