import os
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D


def Get_MSDS_dir(rootdir):
    directory = os.listdir(rootdir)
    directory = [element for element in directory if element != '.DS_Store' and element != 'packed_entropy_data']
    return directory

def load_pickles(path):
    pickle = cPickle.load(open(path))
    return pickle

datafolders = []
root = "/Users/Medina/cellmodeller"
rootdir = root+"/MSD_AAS"
    
filelist = Get_MSDS_dir(rootdir)
nbins = 15
t1 = 400
t2 = 1400

master_msd = np.array([load_pickles(rootdir+'/'+path) for path in filelist])
#shape : (nfiles,r_bins,time)
times = np.arange(t1,t2)
flattened = np.average(master_msd,axis = 0)
data2 = []
for element in master_msd:
    plt.plot(np.log10(times),np.log10(element),linewidth = 0.2)
plt.plot(np.log10(times),np.log10(flattened),"r")
''' 
for data in master_msd:
    for t in data:
        data2.append(t)
pl.hist(data2, bins=np.logspace(np.log10(t1),np.log10(t2), 50))
pl.gca().set_xscale("log")
'''