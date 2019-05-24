import numpy as np
import Conv_curves_lowram as CLR
import Networking as NWK
import cPickle
import matplotlib.pyplot as plt
from RunAndCompress import GetSubDir

def meansqd(array):
    r = array
    n = len(array)
    avg_r = np.average(r)
    r_dev = r-avg_r
    r_sum = np.sum(r_dev**2)
    MSD= r_sum/n
    return MSD
    
root = "/Users/Medina/cellmodeller"
#root = "/media/inmedina/Elements/cellmodeller"
#root = "/home/inmedina/cellmodeller"
datadir = root+"/data"

    
datafolders,datafiles,folders = GetSubDir(datadir)
#print sys.argv[1]
#t1 = int(sys.argv[1])
#t2 = int(sys.argv[2])

#tlist = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
t1 = 150
t2 = 1000
t_tree = 1000

'''
#t2 = 1700 
#t_tree = 1700
path_to_write = root+"/GammaRad"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
'''
i = 0

for simulation in datafiles:

    print 'Loading and running '+ datafolders[i]

    Tree = cPickle.load(open(datafolders[i]+"/Tree_"+str(t_tree)+".pickle","rb"))
    cellstate_0,lin_0 = NWK.loadPickle_lite(simulation,t1)
    bnumber = len(cellstate_0)
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    cellstate_f = CLR.add_radius_angle_area(cellstate_f)
    
    Tree.set_t0_branches(cellstate_0,t1)
    
    MSD_r_array = np.zeros((bnumber,t2-t1))
    MSD_phi_array = np.zeros((bnumber,t2-t1))

    ncells = np.zeros((bnumber,t2-t1))
    
    #need to find new t0 with n cells which i will track the branches
    #ACA ESTA EL PROBLEMA
    bid_array = [bid for bid,branch in Tree.branch.iteritems() if branch.t0 == t1]
    for t in range(t1,t2):
        print 'v----',t
        cellstate,lineage = NWK.loadPickle_lite(simulation,t) 
        cellstate = CLR.add_radius_angle_area(cellstate)
        j = 0
        for bid in bid_array:
        
            branch = Tree.branch[bid]
            r_phi_list = []
            ncell = 0
            
            add = False
            
            branch_cells = [cellstate[node] for node in branch.nodes  if cellstate.has_key(node)]
            r_dist = np.array([cell.r_dist for cell in branch_cells])
            phi = np.array([cell.phi for cell in branch_cells])
            
            for cell_phi in phi:
                if cell_phi >= 7*np.pi/4 or cell_phi <= np.pi/4:
                    add = True
            if add == True:
                for l in range(len(phi)):
                    if phi[l] >= 0 and phi[l] < np.pi: 
                        phi[l] = phi[l]+np.pi
                    elif phi[l] <= 2*np.pi and phi[l] >= np.pi:
                        phi[l] = phi[l]-np.pi
                add = False
                
            ncell = len(branch_cells)
                
            MSD_r = meansqd(r_dist)
            MSD_phi = meansqd(phi)
            
            MSD_r_array[j,t-t1] = MSD_r
            MSD_phi_array[j,t-t1] = MSD_phi
            ncells[j,t-t1] = ncell
            j +=1

    
    