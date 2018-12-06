import os
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import Conv_curves_lowram as CLR
from RunAndCompress import GetSubDir

root = "/Users/Medina/cellmodeller"
#root = "/media/inmedina/Elements/cellmodeller"
#root = "/home/inmedina/cellmodeller"
datadir = root+"/data"
       
    
datafolders,datafiles,folders = GetSubDir(datadir)
    
t1 = 300
t2 = 500
t_tree = 500
#t2 = 1700 
#t_tree = 1700
i = 0
Gamma_hist = []
Nb_hist = []
path_to_write = root+"/histograms"
pw1 = path_to_write+"/Gamma"
pw2 =  path_to_write+"/Nbranch"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
if not os.path.isdir(pw1):
    os.makedirs(pw1)
if not os.path.isdir(pw2):
    os.makedirs(pw2)

for simulation in datafiles:
    Tree = cPickle.load(open(datafolders[i]+"/Tree_"+str(t_tree)+".pickle","rb"))
    
    cellstate_i,lineage_i = CLR.loadPickle_lite(simulation,t1)
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    
    bid_array = np.array([id for id,cell in cellstate_i.iteritems()])
    cell_array = np.array([id for id,cell in cellstate_f.iteritems()])
    
    expected = 1.0/len(cellstate_i)
    N_T = float(len(cellstate_f))
     
    for bid in bid_array:
        N_b = 0.0
        branch_cells = [cellstate_f[id] for id in Tree.branch[bid].nodes if cellstate_f.has_key(id)]
        N_b = float(len(branch_cells))
        Gamma = (N_b/N_T) / expected
        if Gamma != 0:
            Gamma_hist.append(Gamma)
            Nb_hist.append(N_b)
        else:
            print simulation, bid
    i += 1
np.savetxt(pw1+"/Gamma-"+str(t1)+"-"+str(t2)+"-"+str(t_tree)+'.gz',Gamma_hist)
np.savetxt(pw2+"/Nb-"+str(t1)+"-"+str(t2)+"-"+str(t_tree)+'.gz',Nb_hist)
