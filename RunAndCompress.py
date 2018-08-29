import os
import sys
import PickleToGrid as PTG
import cPickle
import AddProteins as App
import numpy as np
import matplotlib.pyplot as plt

def GetSubDir(foldername):
    directory = os.listdir(foldername)
    directory = [element for element in directory if element != '.DS_Store' and element != 'packed_entropy_data']
    filelist = [foldername+'/'+element+'/'+'step-%05d.pickle' for element in directory]
    folderlist = [foldername+'/'+element for element in directory]
    return folderlist,filelist,directory

def obtain_entropy(grid,desired):
    desired_array = np.zeros((grid.gx*grid.gy,grid.nframes,3+len(desired)))
    empty_values = np.array([0.0 for iterator in desired])

    for ix in range(grid.gx):
        for iy in range(grid.gy):
            entropy_of_ensemble = np.array([])
            for it in range(grid.nframes):
                try:
                    desired_values_array = np.array([grid[it,ix,iy].entropy[desired_entropy] for desired_entropy in desired])
                    entropy_of_ensemble = np.append([ix,iy,it],desired_values_array)
                except KeyError:
                    entropy_of_ensemble = np.append([ix,iy,it],empty_values)
                    
                desired_array[ix*grid.gx+iy,it] = entropy_of_ensemble

    return desired_array
            
    
def Compress_grid_entropy(grid,desired):
    comp_info = np.array([grid.nframes,grid.gx,grid.gy,grid.dgx,grid.dgy,grid.center,grid.resize,grid.dt,grid.forwards])
    packed_entropy = obtain_entropy(grid,desired)
    return comp_info,packed_entropy
            


def Calculate_expected(grid):
    
    for it in range(1,grid.nframes-1):
        for ix in range(grid.gx):
            for iy in range(grid.gy):
                Ensemble = grid[it,ix,iy]
                if Ensemble.cell_number > 0.0:  
                    #Calculate the expected value over an Ensemble:
                    #make list of red green etc proteins in an ensemble
                    
                    #logsq_red_prot_array = np.array([Ensemble.cells[id].log_red_squared for (id,cell) in Ensemble.cells.iteritems()])
                    #logsq_green_prot_array = np.array([Ensemble.cells[id].log_green_squared for (id,cell) in Ensemble.cells.iteritems()])
                    #logsq_blue_prot_array = np.array([Ensemble.cells[id].log_blue_squared for (id,cell) in Ensemble.cells.iteritems()])
                    
                    #poissq_red = np.array([Ensemble.cells[id].poisson_red_squared for (id,cell) in Ensemble.cells.iteritems()])
                    
                    Rp_ensemble_array = np.array([cell.red_protein for (id,cell) in Ensemble.cells.iteritems()])
                    vol_ensemble_array = np.array([cell.volume for (id,cell) in Ensemble.cells.iteritems()])
                    concentration_ensemble_array = Rp_ensemble_array/vol_ensemble_array
                    
                    #Ensemble.expect_logsq_red = np.var(logsq_red_prot_array)
                    #Ensemble.expect_logsq_green = np.var(logsq_green_prot_array)
                    #Ensemble.expect_logsq_blue = np.var(logsq_blue_prot_array)
                    #Ensemble.exp_poissq_red = np.var(poissq_red)
                    
                    Ensemble.expected_concentration = np.var(concentration_ensemble_array) 
                    Ensemble.total_Rp = np.sum(Rp_ensemble_array)
                    
                else:
                    Ensemble.expected_concentration = 0.0
                    Ensemble.total_Rp = 0.0
                    #Ensemble.expect_logsq_red = 0.0
                    #Ensemble.expect_logsq_green = 0.0
                    #Ensemble.expect_logsq_blue = 0.0
                    #Ensemble.exp_poissq_red = 0.0
                    
    return grid
                    

def plot_ensemble_all_t(grid,ix,iy):
    max_dist = 0   
     
    array2plot = np.array([grid[it,ix,iy].expected_concentration for it in range(1,grid.nframes-1)])
        
    dist_from_center = np.sqrt(((ix-7.0)**2)+((iy-7.0)**2))
    if dist_from_center > max_dist:
        max_dist = dist_from_center
    color_RGB = [dist_from_center/max_dist,0,0]
    plt.plot(array2plot[::-1],color = color_RGB)
    
def plot_grid(grid):
    for ix in range(grid.gx):
        for iy in range(grid.gy):
            plot_ensemble_all_t(grid,ix,iy)
    
def variance_concentration(cellstates):
    listplot = []

    for t in range(len(cellstates)):
        red_protein_array = np.array([cell.red_protein for (id,cell) in cellstates[t].iteritems()])
        volume_array = np.array([cell.volume for (id,cell) in cellstates[t].iteritems()])
        var_t = np.var(red_protein_array/volume_array)
        print red_protein_array/volume_array
        print var_t
        listplot.append(var_t)
    plt.plot(listplot[::-1],"b")
    
    
'''
def protein_counter(cellstates):
    protein_total_list = []
    for it in range(len(cellstates)):
        red_prot_list =[cell.red_protein for id,cell in cellstates[it].iteritems()]
        print len(cellstates[it]), red_prot_list
        print np.var(red_prot_list)
        protein_total_list.append(np.sum(red_prot_list))
    return protein_total_list[::-1]

def protein_counter_ensemble(grid,ix,iy):
    protein_total_list = []
    for it in range(1,grid.nframes-1):
        
        protein_total_list.append(grid[it,ix,iy].total_Rp)
        
        
    plt.plot(protein_total_list[::-1])
'''        
    
datafolders = []
rootdir = "/Users/Medina/cellmodeller/data"
startframe = 0
nframes = 550
dt = 1 #There's a bit of trouble with this
gridsize = 16
forwards = False
add_proteins = True
#sigma = 0
lambd = 3.0
'''
foldername = sys.argv[1]
startframe = sys.argv[2]
nframes = sys.argv[3]
dt = sys.argv[4] #There's a bit of trouble with this
gridsize = sys.argv[5]
forwards = sys.argv[6]
'''         

datafolders,datafiles,folders = GetSubDir(rootdir)
datapack = []
desired = ['log_red_protein','log_green_protein','log_blue_protein']
i = 0
for simulation in datafiles:
    print 'Loading and running '+ datafolders[i]

    if add_proteins == True:
        cellstates,lineage = App.add_protein_pickles(simulation,startframe,nframes,PTG = forwards,lambd = lambd)
  
    grid,cs = PTG.main(simulation,startframe,nframes,dt,gridsize,forwards = forwards, App = [cellstates,lineage])
    
    grid = Calculate_expected(grid)

    total_var_red = np.array([np.var(np.array([cs[t][id].red_protein for (id,cell) in cs[t].iteritems() ])) for t in range(grid.nframes)])

    print 'Finished, Compressing and writing'
    #comp_info,packed_entropy,packed_variance = Compress_grid(grid,desired)
    
    
    #path_to_write = rootdir+'/packed_entropy_data/'
    #if not os.path.isdir(path_to_write):
    #    os.makedirs(path_to_write)
    #cPickle.dump( packed_entropy, open(path_to_write+folders[i]+".pickle", "wb" ) )
    #i += 1
    
    #print 'Done'


        
        