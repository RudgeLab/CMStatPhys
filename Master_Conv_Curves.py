import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Get_conv_dirs(rootdir):
    directory = os.listdir(rootdir)
    directory = [element for element in directory if element != '.DS_Store' and element != 'packed_entropy_data']
    return directory

def load_arrays(path):
    array = np.loadtxt(path)
    return array

datafolders = []
root = "/Users/Medina/cellmodeller"
rootdir = root+"/convergent_curves"
    
startframe = 0      
    
filelist = Get_conv_dirs(rootdir)
i = 0

nframes = 700
lambd = 1.0
nbins = 25
t1 = 400
t2 =  690


master_conv = np.array([load_arrays(rootdir+'/'+path) for path in filelist])
#shape : (nfiles,delta times,nbins-1)
flattened = np.average(master_conv,axis = 0)
t,r = flattened.shape
r,t = np.arange(r),np.arange(t)
r,t = np.meshgrid(r,t)

'''
ax = Axes3D(plt.gcf())
ax.set_xlabel('r')
ax.set_ylabel('t')
ax.set_zlabel('n_norm')
ax.plot_wireframe(r,t,flattened)
'''

r = flattened.shape[1]
r = np.linspace(0,1,r)
work = flattened[1498]
ababacab = (1-lambd*(r))
plt.plot(ababacab,'k')
plt.plot(work,'r')


