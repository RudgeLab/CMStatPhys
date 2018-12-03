import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Get_MSDS_dir(rootdir):
    directory = os.listdir(rootdir)
    directory = [element for element in directory if element != '.DS_Store' and element != 'packed_entropy_data']
    return directory

def load_arrays(path):
    array = np.loadtxt(path)
    return array

datafolders = []
root = "/Users/Medina/cellmodeller"
rootdir = root+"/MSDS/600"
    
filelist = Get_MSDS_dir(rootdir)
nbins = 15
t1 = 800
tm = 1200
t2 =  1500


master_msd = np.array([load_arrays(rootdir+'/'+path) for path in filelist])
#shape : (nfiles,r_bins,time)
r = np.linspace(0,1,nbins)
times = np.arange(t1,t2)
times_2 = times[times>tm]
flattened = np.average(master_msd,axis = 0)

alphas = []
for variance_r in flattened:
    color = (np.random.rand(1)[0],np.random.rand(1)[0],np.random.rand(1)[0])
    log_time = np.log10(times)
    log_var = np.log10(variance_r)
    
    variances_r_2 = variance_r[times>tm]
    
    if len(variances_r_2)>10:
        polyfit = np.polyfit(np.log10(times_2),np.log10(variances_r_2),1)
        slope = polyfit[0]
        c = polyfit[1]

        #fig, ax = plt.subplots()
        #ax.plot(times, variance_r, 'r', label='MSD')
        #plt.plot(log_time,log_var,"r")
        #plt.plot(np.log10(times_2),slope*np.log10(times_2)+c,"k--")
        alphas.append(slope)
        
fig, ax = plt.subplots()
plt.xlabel("$r/R(t)$")
plt.ylabel("$\\alpha(r)$")
ax.plot(r, alphas, 'r', label='t_r = 400')
plt.grid()
legend = ax.legend(loc='upper left', shadow=True, fontsize='x-large')
