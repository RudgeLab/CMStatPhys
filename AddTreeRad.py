import os
import sys
import numpy as np
import cPickle
import Conv_curves_lowram as CLR
from RunAndCompress import GetSubDir

root = "/Users/Medina/cellmodeller"
#root = "/media/inmedina/Elements/cellmodeller"
#root = "/home/inmedina/cellmodeller"
datadir = root+"/data"
       
    
datafolders,datafiles,folders = GetSubDir(datadir)
#print sys.argv[1]
#t1 = int(sys.argv[1])
#t2 = int(sys.argv[2])

t1 = 0
t2 = 500
t_tree = 500
#t2 = 1700 
#t_tree = 1700
i = 0


for simulation in datafiles:
    path_to_write = datafolders[i]

    print 'Loading and running '+ datafolders[i]

    Tree = cPickle.load(open(datafolders[i]+"/Tree_"+str(t_tree)+".pickle","rb"))
    
    for bid,branch in Tree.branch.iteritems():
        branch.nT = {}
        branch.nB = {}
        branch.r = {}
    
    for t in range(t1,t2+1):
        cellstate,lineage = CLR.loadPickle_lite(simulation,t)
        for id,cell in cellstate.iteritems():
            if Tree.branch[id].t0 == -1:
                print id
                Tree.branch[id].t0 = t
                x_cell = cell.pos[0]
                y_cell = cell.pos[1]
                r_dist = x_cell**2 + y_cell**2
                Tree.branch[id].r0 = np.sqrt(r_dist)            
        for bid,branch in Tree.branch.iteritems():
            Tree.branch[bid].nT[t] = len(cellstate)
            branch_cells = []
            for mid in branch.nodes:
                if cellstate.has_key(mid):
                    branch_cells.append(cellstate[mid])
            branch.nB[t] = len(branch_cells)
                
            
    cPickle.dump(Tree,open(path_to_write+"/Tree_rad_"+str(t_tree)+".pickle","w"))
    i+=1
    print "-----------------------"
