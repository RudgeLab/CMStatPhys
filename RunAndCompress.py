import os
import sys
import PickleToGrid as PTG
import cPickle
import AddProteins as App
import numpy as np
#import matplotlib.pyplot as plt

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
            
    
def Compress_grid(grid,desired):
    comp_info = np.array([grid.nframes,grid.gx,grid.gy,grid.dgx,grid.dgy,grid.center,grid.resize,grid.dt,grid.forwards])
    packed_entropy = obtain_entropy(grid,desired)
    return comp_info,packed_entropy
            


def Calculate_expected_squared(grid):
    
    for it in range(1,grid.nframes-1):
        for ix in range(grid.gx):
            for iy in range(grid.gy):
                Ensemble = grid[it,ix,iy]
                if Ensemble.cell_number > 0.0:  
                    #Calculate the expected value over an Ensemble:
                    for id,cell in Ensemble.cells.iteritems():
                        
                        cell.log_red_squared =  cell.log_red_protein
                        cell.log_green_squared = cell.log_green_protein
                        cell.log_blue_squared = cell.log_blue_protein
                        
                    #make list of red green etc proteins in an ensemble
                    
                    logsq_red_prot_array = np.sort(np.array([Ensemble.cells[id].log_red_squared for (id,cell) in Ensemble.cells.iteritems()]))
                    logsq_green_prot_array = np.sort(np.array([Ensemble.cells[id].log_green_squared for (id,cell) in Ensemble.cells.iteritems()]))
                    logsq_blue_prot_array = np.sort(np.array([Ensemble.cells[id].log_blue_squared for (id,cell) in Ensemble.cells.iteritems()]))
                    
                    Ensemble.expect_logsq_red = np.var(logsq_red_prot_array)
                    Ensemble.expect_logsq_green = np.var(logsq_green_prot_array)
                    Ensemble.expect_logsq_blue = np.var(logsq_blue_prot_array)
                else:
                    Ensemble.expect_logsq_red = 0.0
                    Ensemble.expect_logsq_green = 0.0
                    Ensemble.expect_logsq_blue = 0.0
    return grid
                    

def plot_ensemble_all_t(grid,ix,iy):
    listplot = []
    for it in range(nframes):
        listplot = listplot.append(grid[it,ix,iy].expect_logsq_red)
    plt.plot(listplot)
    
def plot_grid(grid):
    for ix in range(grid.gx):
        for iy in range(grid.gy):
            plot_ensemble_all_t(grid,ix,iy)
    

    
    
datafolders = []
rootdir = "/Users/Medina/cellmodeller/data"
startframe = 0
nframes = 130
dt = 1 #There's a bit of trouble with this
gridsize = 16
forwards = False
add_proteins = True
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
        cellstates,lineage = App.add_protein_pickles(simulation,startframe,nframes,PTG = forwards)
  
    grid,cs = PTG.main(simulation,startframe,nframes,dt,gridsize,forwards = forwards, App = [cellstates,lineage])
    
    
    
    grid = Calculate_expected_squared(grid)

    total_var_red = np.array([np.var(np.array([cs[t][id].log_red_protein for (id,cell) in cs[t].iteritems() ])) for t in range(grid.nframes)])

    print 'Finished, Compressing and writing'
    comp_info,packed_entropy = Compress_grid(grid,desired)
    
    
    path_to_write = rootdir+'/packed_entropy_data/'
    if not os.path.isdir(path_to_write):
        os.makedirs(path_to_write)
    cPickle.dump( packed_entropy, open(path_to_write+folders[i]+".pickle", "wb" ) )
    i += 1
    
    print 'Done'


        
        