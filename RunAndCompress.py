import os
import sys
import PickleToGrid as PTG
import cPickle
import numpy as np

def GetSubDir(foldername):
    directory = os.listdir(foldername)
    directory = [element for element in directory if element != '.DS_Store']
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
                
        
    

    
    
datafolders = []
rootdir = "/Users/Medina/cellmodeller/data"
startframe = 0
nframes = 180
dt = 1 #There's a bit of trouble with this
gridsize = 16
forwards = False

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
    grid,cs = PTG.main(simulation,startframe,nframes,dt,gridsize,forwards = forwards)
    
    print 'Finished, Compressing and writing'
    comp_info,packed_entropy = Compress_grid(grid,desired)
    
    path_to_write = rootdir+'/packed_entropy_data/'
    if not os.path.isdir(path_to_write):
        os.makedirs(path_to_write)
    cPickle.dump( packed_entropy, open(path_to_write+folders[i]+".pickle", "wb" ) )
    i += 1
    
    print 'Done'


        
        