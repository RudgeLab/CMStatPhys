import Networking as NWK
import BranchDistributionnew as BDN
import Conv_curves_lowram as CLR
import numpy as np
import os
from RunAndCompress import GetSubDir


root = "/Users/Medina/cellmodeller"
#root = "/home/inmedina/cellmodeller"
rootdir = root+"/data"
       
    
datafolders,datafiles,folders = GetSubDir(rootdir)



t1 = 400
tm = 600
t2 = 1700
nbins = 15
i = 0
path_to_write = root+"/MSDS"

if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
    
for simulation in datafiles:
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    print 'Loading and running '+ datafolders[i]
    Willow = NWK.Create_Tree(lineage_f)
    Willow.create_all_branches()
    print "obtaining MSD and writing"
    r,MSDS,times,alphas = BDN.run_script(simulation,t1,t2,tm,Willow,nbins)
    
    np.savetxt(path_to_write+'/'+folders[i]+str(nbins)+'.gz',MSDS)
    i+=1
    print "-----------------------"
    
    