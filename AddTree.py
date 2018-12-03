import Networking as NWK
import Conv_curves_lowram as CLR
import numpy as np
import os
import cPickle
from RunAndCompress import GetSubDir


root = "/Users/Medina/cellmodeller"
#root = "/home/inmedina/cellmodeller"
datadir = root+"/data"
       
    
datafolders,datafiles,folders = GetSubDir(datadir)

t2 = 1700
i = 0
for simulation in datafiles:
    path_to_write = datafolders[i]
    
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    
    print 'Loading and running '+ datafolders[i]
    Willow = NWK.Create_Tree(lineage_f)
    Willow.create_all_branches()
    
    print "obtaining MSD and writing"    
    cPickle.dump(Willow,open(path_to_write+"/Tree_"+str(t2)+".pickle","wb"))
    i+=1
    print "-----------------------"
