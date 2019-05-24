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
t1 = 0
t2 = 700
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
    path_to_write = datafolders[i]

    print 'Loading and running '+ datafolders[i]

    Tree = cPickle.load(open(datafolders[i]+"/Tree_"+str(t_tree)+".pickle","rb"))
    cellstate_0,lin_0 = NWK.loadPickle_lite(simulation,t1)
    bnumber = len(cellstate_0)
    cellstate_f,lineage_f = CLR.loadPickle_lite(simulation,t2)
    cellstate_f = CLR.add_radius_angle_area(cellstate_f)
    
    #Tree.set_t0_branches(cellstate_0,t1)
    sim_cells_phi = {}
    sim_cells_t0 = {}
    sim_cells_n ={} 
    sim_cells_r = {}
    sim_cells_t = {}
    for t in range(t1,t2):
        print 'v----',t
        cellstate,lineage = NWK.loadPickle_lite(simulation,t) 
        cellstate = CLR.add_radius_angle_area(cellstate)
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
           
            
    #determinando l_t = r_t*(theta(t)-theta(t_0))
    sim_cells_rtheta = {}   
    init_times = {}
    plt.figure()
    for id,list in sim_cells_phi.iteritems():
        array = np.array(list)
        array_init = array - array[0]      
        rdev = sim_cells_r[id]*array_init
        sim_cells_rtheta[id] = rdev
        
        init_times[id] = np.array(sim_cells_t[id])-sim_cells_t[id][0]
        plt.plot(init_times[id],sim_cells_rtheta[id],"bo",markersize = 0.2)
     
        
    #binning in time
    bin_dic = {}
    for id,time_array in init_times.iteritems():
        for t in time_array:
            if t not in bin_dic.keys():
                bin_dic[t] = []
            bin_dic[t].append(sim_cells_rtheta[id][t])
            
    #MSD bins:
    MSD_master = []
    plt.figure()
    for time,bin_array in bin_dic.iteritems():
        MSD = meansqd(bin_array)
        MSD_master.append(MSD)
        plt.plot(np.log10(time),np.log10(MSD),"bo",markersize = 0.2)
        
    cPickle.dump(MSD_master,open(path_to_write+"/MSD_AA_"+str(t1)+"-"+str(t2)+".pickle","w"))
    i+=1
    print "-----------------------"
        
        
        
'''
    all_msd = {}
    for id,list in sim_cells_phi.iteritems():
        array = np.array(list)
        array_init = array - array[0]
        MSD = meansqdr(array_init,sim_cells_r[id])
        all_msd[id] = MSD
    for id,msd in all_msd.iteritems():
        
        plt.plot(sim_cells_t[id],np.log10(msd),"ro",markersize=0.2) #correr esto sin modificar
        '''
        
        
        
'''    
    MSD_r_array = np.zeros((bnumber,t2-t1))
    MSD_phi_array = np.zeros((bnumber,t2-t1))

    ncells = np.zeros((bnumber,t2-t1))
    
'''
    
    