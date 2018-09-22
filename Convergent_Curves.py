import PlayRoom as PR
import AddProteins as App
from RunAndCompress import GetSubDir
import os
import numpy as np

datafolders = []
root = "/Users/Medina/cellmodeller"
rootdir = root+"/data"
    
startframe = 0      
    
datafolders,datafiles,folders = GetSubDir(rootdir)
i = 0

nframes = 700
lambd = 1.0
nbins = 20
t1 = 10
t2 =  690
repetitions = 3
path_to_write = root+"/convergent_curves"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
for simulation in datafiles:
    for repet in range(1,repetitions+1):
        print "Repetition: ", repet
        print 'Loading and running '+ datafolders[i]
        cellstates = App.add_protein_pickles(simulation,startframe,nframes,lambd = lambd)
        cellstates_reordered = [cellstates[t]['cellStates'] for t in range(nframes)]
    
        cellstate = App.add_radius(cellstates_reordered) #get r for each cell
    
        n_norm, r_bins = PR.obtain_convergent_curves(cellstate,t1,t2,nbins)
    
        np.savetxt(path_to_write+'/'+folders[i]+'-r-'+str(repet)+'.gz',n_norm)
    i+=1
    
