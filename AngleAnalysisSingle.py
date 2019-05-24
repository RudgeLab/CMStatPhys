import numpy as np
import Conv_curves_lowram as CLR
import Networking as NWK
import cPickle
import matplotlib.pyplot as plt
import os
from RunAndCompress import GetSubDir

def meansqd(array):
    r = array
    n = len(array)
    avg_r = np.average(r)
    r_dev = r-avg_r
    r_sum = np.sum(r_dev**2)
    MSD= r_sum/n
    return MSD

def meansqdr(array,rdist):
    r = array
    n = len(array)
    avg_r = np.average(r)
    r_dev = r-avg_r
    r_sum = np.sum((rdist*r_dev)**2)
    MSD= np.sqrt(r_sum/n)
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
t1 = 255
t2 = 700
t_tree = 500

'''
#t2 = 1700 
#t_tree = 1700
path_to_write = root+"/GammaRad"
if not os.path.isdir(path_to_write):
    os.makedirs(path_to_write)
'''
i = 0

for simulation in datafiles:
    path_to_write = datafolders[i]
    print 'Loading and running '+ datafolders[i]
    cellstate_0,lin_0 = NWK.loadPickle_lite(simulation,t1)
    bnumber = len(cellstate_0)
    
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    cellstate_f = CLR.add_radius_angle_area(cellstate_f)
    
    cell_select_list = [id for id,cell in cellstate_0.iteritems()]
    print cell_select_list
        
    
    #Tree.set_t0_branches(cellstate_0,t1)
    sim_cells_phi = {}
    sim_cells_t0 = {}
    sim_cells_n ={} 
    sim_cells_r = {}
    sim_cells_t = {}
    
    sim_cells_id = {}
    for t in range(t1,t2):
        print 'v----',t
        cellstate,lineage = NWK.loadPickle_lite(simulation,t) 
        sim_cells_id[t] = [id for id,cell in cellstate.iteritems()]
        cellstate = CLR.add_radius_angle_area(cellstate)
        inv_map = {}
        for did,pid in lineage.iteritems():
            inv_map[pid] = inv_map.get(pid, [])
            inv_map[pid].append(did)
            
        for id,cell in cellstate.iteritems():
            if id not in sim_cells_phi.keys():
                sim_cells_phi[id] = []
                sim_cells_t0[id]= t
                sim_cells_t[id] = []
                sim_cells_n[id] = len(cellstate)
                sim_cells_r[id] = []
            sim_cells_phi[id].append(cell.phi)
            sim_cells_r[id].append(cell.r_dist)
            sim_cells_t[id].append(t)

            #plt.plot(t,np.pi-cell.phi,"ro",markersize=0.4)
    new_cells_phi = {}
    new_cells_r = {}
    for cell_select in cell_select_list:
        new_id = cell_select
        for t in range(t1,t2):

            if cell_select in sim_cells_id[t] and cell_select not in new_cells_phi.keys():
                new_cells_phi[cell_select] = sim_cells_phi[cell_select]
                new_cells_r[cell_select] = sim_cells_r[cell_select]
                new_id = cell_select
                
            elif new_id not in sim_cells_id[t]:
                choices = inv_map[new_id]
                new_id = choices[np.random.randint(0,2)]
                if new_id in sim_cells_id[t]:
                    new_cells_phi[cell_select] += sim_cells_phi[new_id]
                    new_cells_r[cell_select] += sim_cells_r[new_id]

    new_cells_phi_temp = dict(new_cells_phi)
    
    for id,listt in new_cells_phi_temp.iteritems():
        if any(element > 5.8 for element in listt) and any(element < 0.5 for element in listt):
            new_cells_phi.pop(id,None)
            new_cells_r.pop(id,None)
            print "happens"
    
         
    #determinando l_t = r_t*(theta(t)-theta(t_0))
    sim_cells_rtheta = {}   
    for id,list in new_cells_phi.iteritems():
        array = np.array(list)
        array_init = array - array[0]      
        rdev = new_cells_r[id]*array_init
        sim_cells_rtheta[id] = rdev
        
    #plt.figure()
    #for id,rthings in sim_cells_rtheta.iteritems():
        #plt.plot(rthings)
    #binning in time
    
    bin_dic = {}
    time_array = np.arange(0,len(sim_cells_rtheta[id]))
    for id,rtheta in sim_cells_rtheta.iteritems():
        for t in time_array:
            if t not in bin_dic.keys():
                bin_dic[t] = []
            bin_dic[t].append(rtheta[t])
            
    #MSD bins:
    plt.figure()
    MSD_master = []
    for time,bin_array in bin_dic.iteritems():
        MSD = meansqd(bin_array)
        MSD_master.append(MSD)
        plt.plot(np.log10(time),np.log10(MSD),"bo",markersize = 0.2)
    
    cPickle.dump(MSD_master,open(path_to_write+"/MSD_AAS_"+str(t1)+"-"+str(t2)+".pickle","w"))
    i+=1
    print "-----------------------"

        
        