import os
import numpy as np
import DivisionDistribution as DD
from PickleToGrid import loadPickle_pro

def GetSubDir(foldername):
    directory = os.listdir(foldername)
    directory = [element for element in directory if element != '.DS_Store' and element != 'packed_entropy_data']
    filelist = [foldername+'/'+element+'/'+'step-%05d.pickle' for element in directory]
    folderlist = [foldername+'/'+element for element in directory]
    return folderlist,filelist,directory


def add_proteins(cellstate,lineage,sigma = 0, lambd = 0): #cellstate[t],cellstate[t+1],lineage[t+1]
    '''POSSIBLE PROBLEM: cells update their divideTag in CellModeller simulations
    so its not a good checker for division, ask Tim what other possible checkers exist, 
    hopefully you don't have to mess with the simulations
    '''
    
    #------------- dividing module
    for t in range(len(cellstate)-1):
        for id,parent in cellstate[t].iteritems():
            if t == 0:
                
                #init RGB_log
                #parent.log_red_protein = 1.0
                #parent.log_green_protein = 1.0
                #parent.log_blue_protein = 1.0
                #init RGB_poisson
                #parent.red_poisson = 0.0
                #parent.green_poisson = 0.0
                #parent.blue_poisson = 0.0
                
                #init binomial+poisson
                parent.red_protein = 0.0
                
            #print cellstate[t]
            #print lineage[t]
            #print lineage[t+1]
            try:
                #print "Testing", parent
                cellstate[t+1][id] #divides?
                
                #print parent.divideFlag #Doesn't work in CM

                if parent.divideFlag == False:
                #update not divided
                    #print cellstate[t+1]
                    cellstate[t+1][id].red_protein = parent.red_protein + np.random.poisson(lambd)
                    '''
                    cellstate[t+1][id].log_red_protein = cellstate[t][id].log_red_protein
                    cellstate[t+1][id].log_green_protein = cellstate[t][id].log_green_protein
                    cellstate[t+1][id].log_blue_protein = cellstate[t][id].log_blue_protein
                    
                    cellstate[t+1][id].red_poisson = cellstate[t][id].red_poisson
                    cellstate[t+1][id].green_poisson = cellstate[t][id].green_poisson
                    cellstate[t+1][id].blue_poisson = cellstate[t][id].blue_poisson
                    '''
            #division
            except KeyError:
                #print "KE"
                #print parent.divideFlag
                                    
                dids = [key for key, value in lineage[t+1].iteritems() if value == id]
                did1 = dids[0]
                did2 = dids[1]
                #print 'dids: ', did1, did2
                #comment unwanted modules:
                #-------Normal Log
                #cellstate[t][id],cellstate[t+1][did1],cellstate[t+1][did2] = DD.divide_RGB_log(parent,cellstate[t+1][did1],cellstate[t+1][did2],sigma)
                #------Poisson distribution
                #cellstate[t][id],cellstate[t+1][did1],cellstate[t+1][did2] = DD.divide_RGB_Poisson(parent,cellstate[t+1][did1],cellstate[t+1][did2],lambd)
                #------Binomial division
                cellstate[t][id],cellstate[t+1][did1],cellstate[t+1][did2] = DD.divide_binomial(parent,cellstate[t+1][did1],cellstate[t+1][did2])
            #print '+++++++++++next cell'
    return cellstate

def add_protein_pickles(fname,startframe,nframes,forwards = True, PTG = None,sigma = 0.5, lambd = 5):
    simplepickle = []
    
    cellstate, lineage = loadPickle_pro(fname,startframe,nframes,1,True)
    
    cellstate = add_proteins(cellstate,lineage, lambd=lambd)
    
    for t in range(nframes):
        auxdic = {}
        auxdic['cellStates'] = cellstate[t]
        auxdic['lineage'] = lineage[t]
        simplepickle.append(auxdic)
        
    if PTG == True:
        return cellstate,lineage
    
    if PTG == False:
        return cellstate[::-1],lineage[::-1]
    
    else:
        return simplepickle
    
def add_radius(cellstate):
    for it in range(len(cellstate)):
        for id,cell in cellstate[it].iteritems():
            x_cell,y_cell,z_cell = cell.pos[0],cell.pos[1],cell.pos[2]
            cellstate[it][id].r_dist =  np.sqrt((x_cell**2)+(y_cell**2)+(z_cell**2))
    return cellstate
            
    
'''   
datafolders = []
rootdir = "/Users/Medina/cellmodeller/data"
startframe = 0
nframes = 180
dt = 1 #There's a bit of trouble with this
gridsize = 16
forwards = False
full_folders,data_names_files,folders = GetSubDir(rootdir)

i = 0
for simulation in datafiles:
    newpickle = add_proteins(simulation,startframe,nframes,forwards = forwards)
    
    path_to_write = rootdir+'/rewritten_pickles/'
    
    
    if not os.path.isdir(path_to_write):
        os.makedirs(path_to_write)
    
    for t in range(len(newpickle)):
        #for now, after don't use writing
        new_path_to_write = path_to_write+folders[i]
        if not os.path.isdir(new_path_to_write):
            os.makedirs(new_path_to_write)
        print "Writing: "+new_path_to_write+"/step-%05d"%t+".pickle"
        cPickle.dump( newpickle[t], open(new_path_to_write+"/step-%05d"%t+".pickle", "wb" ) )
    i+=1
    '''
    