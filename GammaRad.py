import os
import sys
import numpy as np
import cPickle
import Conv_curves_lowram as CLR
import matplotlib.pyplot as plt
from RunAndCompress import GetSubDir

def get_R_max_t(cells_t):     #get max R

    R_max_t = 0.0
    for id,cell in cells_t.iteritems():
        r_dist = get_cell_r(cell)
        if r_dist > R_max_t:
            R_max_t = r_dist
    return R_max_t

def get_cell_r(cell):
    x_cell = cell.pos[0]
    y_cell = cell.pos[1]
    r_dist = np.sqrt((x_cell**2) + (y_cell**2))
    return r_dist
    

root = "/Users/Medina/cellmodeller"
#root = "/media/inmedina/Elements/cellmodeller"
#root = "/home/inmedina/cellmodeller"
datadir = root+"/data"
       
    
datafolders,datafiles,folders = GetSubDir(datadir)
#print sys.argv[1]
#t1 = int(sys.argv[1])
#t2 = int(sys.argv[2])

#tlist = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
tlist = [100,200,300]
t2 = 500
t_tree = 500

#t2 = 1700 
#t_tree = 1700
i = 0
path_to_write = root+"/GammaRad"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)

for simulation in datafiles:

    print 'Loading and running '+ datafolders[i]

    Tree = cPickle.load(open(datafolders[i]+"/Tree_"+str(t_tree)+".pickle","rb"))
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    N_T = float(len(cellstate_f))
    for t in tlist:
        print t
        cellstate,lineage = CLR.loadPickle_lite(simulation,t)
        
        R_max_t = get_R_max_t(cellstate)
        expected = 1.0/len(cellstate)
        
        Gamma_list = []
        R_list = []
        for id,cell in cellstate.iteritems():
            #R:
            R = get_cell_r(cell)
            R_list.append(R)
            #Gamma:
            branch = Tree.branch[id]
            branch_cells = [cell for nid,cell in cellstate_f.iteritems() if nid in branch.nodes]
            N_b = float(len(branch_cells))
            N_norm = N_b/N_T
            Gamma = N_norm / expected
            
            Gamma_list.append(Gamma)
            Gamma_array = np.array(Gamma_list)
        
        R_array = np.array(R_list)
        R_norm = R_array/R_max_t
        RminR = R_max_t - np.array(R_list)
        
        #sort:
        
        indexes = range(len(RminR))
        indexes.sort(key=RminR.__getitem__)
        sorted_R = map(RminR.__getitem__, indexes)
        sorted_Gamma = map(Gamma_array.__getitem__, indexes)
        final = np.zeros((2,len(sorted_R)))
        final[0] = sorted_Gamma
        final[1] = sorted_R
        
        np.savetxt(path_to_write+"/GammaRad-"+str(t)+"-"+str(t2)+"-"+str(t_tree)+'.gz',final)
    i+=1
    print "-----------------------"
